#include "Stepper.h"
#include "ModelA.h"
#include "O4AlgebraHelper.h"
#include <algorithm>
#include <cstdio>

///////////////////////////////////////////////////////////////////////////

IdealPV2::IdealPV2(ModelA &in, const bool &accept_reject_in)
    : model(&in), da(model->domain), data(model->data),
      accept_reject(accept_reject_in), monitor("Statistics of IdealPV2") {

  // Store for later use
  DMDAGetCorners(da, &xstart, &ystart, &zstart, &xdimension, &ydimension,
                 &zdimension);
  if (accept_reject) {
    VecDuplicate(model->solution, &previoussolution);
  }
  oldEnergy = 0.;
  newEnergy = 0.;
}

bool IdealPV2::step(const double &dt) {
  bool success = true;

  // Do some setup for the accept reject procedure
  if (accept_reject) {
    PetscCall(VecCopy(model->solution, previoussolution));
  }
  oldEnergy = computeEnergy(dt);

  // Take the proposal
  success = step_no_reject(dt);

  // Do a metropolis accept reject for the ideal step.  Even if we are not
  // doing the accept reject we keep track of the statistics.
  newEnergy = computeEnergy(dt);
  PetscScalar deltaE = newEnergy - oldEnergy;

  PetscBool reject(PETSC_FALSE);
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if (myRank == 0) {
    if (deltaE > 0) {
      PetscScalar prob = ModelARndm->uniform();
      if (prob > exp(-deltaE)) {
        reject = PETSC_TRUE;
        monitor.increment_up_no(deltaE);
      } else {
        monitor.increment_up_yes(deltaE);
      }
    } else {
      monitor.increment_down(deltaE);
    }
  }

  // If we doing the accept reject, we need to acctualy reject
  if (accept_reject) {
    MPI_Bcast(&reject, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    if (reject) {
      PetscCall(VecCopy(previoussolution, model->solution));
    }
  }
  return success;
}

bool IdealPV2::step_no_reject(const double &dt) {

  // drifts the solution by dt / 2.0
  G_node ***phinew;
  PetscCall(DMDAVecGetArray(da, model->solution, &phinew));
  rotatePhi(phinew, dt / 2.0);
  PetscCall(DMDAVecRestoreArray(da, model->solution, &phinew));

  // Get a local vector with ghost cells
  Vec localU;
  PetscCall(DMGetLocalVector(da, &localU));

  // Fill in the ghost celss with mpicalls
  PetscCall(DMGlobalToLocalBegin(da, model->solution, INSERT_VALUES, localU));
  PetscCall(DMGlobalToLocalEnd(da, model->solution, INSERT_VALUES, localU));

  // Get Access to arrays with drifted solution
  G_node ***phi;
  PetscCall(DMDAVecGetArrayRead(da, localU, &phi));
  PetscCall(DMDAVecGetArray(da, model->solution, &phinew));

  const auto &coeff = data.acoefficients;
  PetscReal H[4] = {coeff.H, 0., 0., 0.};
  if (data.ahandler.superfluidmode) {
    H[0] *= coeff.sigmabyf(data.atime.t());
  }
  const PetscReal axx = pow(1. / data.hX(), 2);
  const PetscReal ayy = pow(1. / data.hY(), 2);
  const PetscReal azz = pow(1. / data.hZ(), 2);

  PetscScalar advxx, advyy, advzz;
  PetscInt s1, s2, epsilon;

  // Loop over central elements
  for (PetscInt k = zstart; k < zstart + zdimension; k++) {
    for (PetscInt j = ystart; j < ystart + ydimension; j++) {
      for (PetscInt i = xstart; i < xstart + xdimension; i++) {
        // First evolve the momenta nab

        G_node &centralPhi = phi[k][j][i];
        G_node &phixplus = phi[k][j][i + 1];
        G_node &phixminus = phi[k][j][i - 1];
        G_node &phiyplus = phi[k][j + 1][i];
        G_node &phiyminus = phi[k][j - 1][i];
        G_node &phizplus = phi[k + 1][j][i];
        G_node &phizminus = phi[k - 1][j][i];

        for (PetscInt l = 0; l < ModelAData::NA; l++) {
          advxx = (-phixplus.f[0] * centralPhi.f[l + 1] +
                   centralPhi.f[0] * phixplus.f[l + 1] +
                   centralPhi.f[0] * phixminus.f[l + 1] -
                   phixminus.f[0] * centralPhi.f[l + 1]) *
                  axx;

          advyy = (-phiyplus.f[0] * centralPhi.f[l + 1] +
                   centralPhi.f[0] * phiyplus.f[l + 1] +
                   centralPhi.f[0] * phiyminus.f[l + 1] -
                   phiyminus.f[0] * centralPhi.f[l + 1]) *
                  ayy;

          advzz = (-phizplus.f[0] * centralPhi.f[l + 1] +
                   centralPhi.f[0] * phizplus.f[l + 1] +
                   centralPhi.f[0] * phizminus.f[l + 1] -
                   phizminus.f[0] * centralPhi.f[l + 1]) *
                  azz;

          phinew[k][j][i].A[l] +=
              dt * (advxx + advyy + advzz - H[0] * centralPhi.f[l + 1]);
        }

        for (PetscInt s = 0; s < ModelAData::NV; s++) {

          s1 = (s + 1) % 3;
          s2 = (s + 2) % 3;
          epsilon = ((PetscScalar)(s - s1) * (s1 - s2) * (s2 - s)) / 2.;
          // advection term with epsilon
          advxx = epsilon *
                  (phixplus.f[1 + s2] * centralPhi.f[1 + s1] -
                   centralPhi.f[1 + s2] * phixminus.f[1 + s1] -
                   phixplus.f[1 + s1] * centralPhi.f[1 + s2] +
                   centralPhi.f[1 + s1] * phixminus.f[1 + s2]) *
                  axx;

          advyy = epsilon *
                  (phiyplus.f[1 + s2] * centralPhi.f[1 + s1] -
                   centralPhi.f[1 + s2] * phiyminus.f[1 + s1] -
                   phiyplus.f[1 + s1] * centralPhi.f[1 + s2] +
                   centralPhi.f[1 + s1] * phiyminus.f[1 + s2]) *
                  ayy;

          advzz = epsilon *
                  (phizplus.f[1 + s2] * centralPhi.f[1 + s1] -
                   centralPhi.f[1 + s2] * phizminus.f[1 + s1] -
                   phizplus.f[1 + s1] * centralPhi.f[1 + s2] +
                   centralPhi.f[1 + s1] * phizminus.f[1 + s2]) *
                  azz;

          phinew[k][j][i].V[s] += dt * (advxx + advyy + advzz);
        }
      }
    }
  }

  // drifts dt / 2.0
  rotatePhi(phinew, dt / 2.0);

  PetscCall(DMDAVecRestoreArray(da, model->solution, &phinew));
  PetscCall(DMDAVecRestoreArrayRead(da, localU, &phi));
  PetscCall(DMRestoreLocalVector(da, &localU));

  return true;
}

void IdealPV2::rotatePhi(G_node ***phinew, double dt) {

  const auto &coeff = data.acoefficients;
  for (PetscInt k = zstart; k < zstart + zdimension; k++) {
    for (PetscInt j = ystart; j < ystart + ydimension; j++) {
      for (PetscInt i = xstart; i < xstart + xdimension; i++) {

        // here we have to convert the charge from the chemical potential
        PetscScalar axialmu[ModelAData::NA], vectormu[ModelAData::NV];

        for (PetscInt s = 0; s < ModelAData::NV; s++) {

          vectormu[s] = -phinew[k][j][i].V[s] / coeff.chi * dt;
        }

        for (PetscInt s = 0; s < ModelAData::NA; s++) {

          axialmu[s] = -phinew[k][j][i].A[s] / coeff.chi * dt;
        }

        O4AlgebraHelper::O4Rotation(vectormu, axialmu, phinew[k][j][i].f);
      }
    }
  }
}

PetscScalar IdealPV2::computeEnergy(double dt) {
  DM da = model->domain;
  // Get a local vector with ghost cells
  Vec localUNew;
  DMGetLocalVector(da, &localUNew);

  // Fill in the ghost celss with mpicalls

  DMGlobalToLocalBegin(da, model->solution, INSERT_VALUES, localUNew);
  DMGlobalToLocalEnd(da, model->solution, INSERT_VALUES, localUNew);

  const auto &data = model->data;
  const auto &coeff = data.acoefficients;

  G_node ***phiNew;
  DMDAVecGetArrayRead(da, localUNew, &phiNew);

  const PetscReal H[4] = {coeff.H, 0., 0., 0.};

  PetscInt xstart, ystart, zstart, xdimension, ydimension, zdimension;
  DMDAGetCorners(da, &xstart, &ystart, &zstart, &xdimension, &ydimension,
                 &zdimension);

  // Loop over central elements
  PetscScalar phimid = 0, grad2 = 0, nab2 = 0, hEn = 0;
  PetscScalar localEnergy = 0;
  // std::array<PetscScalar,3> localEnergyArr {0,0,0};

  for (PetscInt k = zstart; k < zstart + zdimension; k++) {
    for (PetscInt j = ystart; j < ystart + ydimension; j++) {
      for (PetscInt i = xstart; i < xstart + xdimension; i++) {
        for (PetscInt s = 0; s < ModelAData::Nphi; ++s) {

          phimid = phiNew[k][j][i].f[s];

          grad2 += pow(phiNew[k + 1][j][i].f[s] - phimid, 2);
          grad2 += pow(phiNew[k][j + 1][i].f[s] - phimid, 2);
          grad2 += pow(phiNew[k][j][i + 1].f[s] - phimid, 2);

          hEn += phimid * H[s];
        }
        for (PetscInt s = 0; s < ModelAData::NV; s++) {
          nab2 += pow(phiNew[k][j][i].V[s], 2);
        }

        for (PetscInt s = 0; s < ModelAData::NA; s++) {
          nab2 += pow(phiNew[k][j][i].A[s], 2);
        }
      }
    }
  }

  localEnergy += 0.5 / coeff.chi * nab2;
  localEnergy += 0.5 * grad2;
  localEnergy -= hEn;

  PetscScalar energy = 0;
  MPI_Reduce(&localEnergy, &energy, 1, MPIU_SCALAR, MPI_SUM, 0,
             PETSC_COMM_WORLD);

  DMDAVecRestoreArrayRead(da, localUNew, &phiNew);
  DMRestoreLocalVector(da, &localUNew);

  return energy;
}

void IdealPV2::finalize() {
  if (accept_reject) {
    VecDestroy(&previoussolution);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if (myRank == 0) {
      monitor.print(stdout);
    }
  }
}
/////////////////////////////////////////////////////////////////////////

EulerLangevinHB::EulerLangevinHB(ModelA &in)
    : model(&in), monitor("Phi HB Steps") {
  DMCreateLocalVector(model->domain, &phi_local);
}

bool EulerLangevinHB::step(const double &dt) {

  // Get the ranges
  PetscInt ixs, iys, izs, nx, ny, nz;
  PetscCall(DMDAGetCorners(model->domain, &ixs, &iys, &izs, &nx, &ny, &nz));

  G_node heff = {};
  G_node phi_o = {};
  G_node phi_n = {};

  const auto &ahandler = model->data.ahandler;
  const auto &atime = model->data.atime;
  const auto &coeff = model->data.acoefficients;
  const double &Lambda2 = 0.5 * (2. * 3. + model->data.mass()); // mass term
  const double &lambda = coeff.lambda;

  double H = coeff.H;
  double rdtg = sqrt(2. * dt * coeff.gamma);

  // In superfluid mode the coefficients are modified.  rdtg,  which is a
  // measure of the field diffusion,  is divided by f2. Similarly the field, H,
  // is multiplied  by a factor $\bar \sigma/f$.
  double f = 0.;
  if (ahandler.superfluidmode) {
    double sigmabyf = coeff.sigmabyf(atime.t());
    H *= sigmabyf;
  }

  PetscLogEvent communication, random, loop;
  PetscLogEventRegister("Communication", 0, &communication);
  PetscLogEventRegister("Random", 0, &random);
  PetscLogEventRegister("Loop", 0, &loop);

  // Checkerboard order ieo = even and odd sites
  for (PetscInt ieo = 0; ieo < 2; ieo++) {
    // take the global vector U and distribute to the local vector localU
    PetscLogEventBegin(communication, 0, 0, 0, 0);
    PetscCall(DMGlobalToLocalBegin(model->domain, model->solution,
                                   INSERT_VALUES, phi_local));
    PetscCall(DMGlobalToLocalEnd(model->domain, model->solution, INSERT_VALUES,
                                 phi_local));
    PetscLogEventEnd(communication, 0, 0, 0, 0);
    // Get pointer to local array
    G_node ***phi;
    PetscCall(DMDAVecGetArray(model->domain, phi_local, &phi));
    // Get pointer global array
    G_node ***phinew;
    PetscCall(DMDAVecGetArray(model->domain, model->solution, &phinew));

    PetscLogEventBegin(loop, 0, 0, 0, 0);
    for (PetscInt k = izs; k < izs + nz; k++) {
      for (PetscInt j = iys; j < iys + ny; j++) {
        for (PetscInt i = ixs; i < ixs + nx; i++) {
          if ((k + j + i) % 2 != ieo) {
            continue;
          }

          PetscScalar hphi_o = 0.; // old and new h*phi
          PetscScalar hphi_n = 0.; // old and new h*phi
          PetscScalar s_o = 0.;    // old sum of phi**2
          PetscScalar s_n = 0.;    // new sum of phi**2

          if (ahandler.superfluidmode) {
            // This is the update for superfluid mode
            for (PetscInt l = 0; l < ModelAData::Nphi; l++) {
              phi_o.f[l] = phi[k][j][i].f[l]; // Fill up the old values
              phi_n.f[l] = phi[k][j][i].f[l]; // Will hold the new values
            }

            // Fill up random angles and do the rotation
            PetscScalar V[3], A[3];
            for (PetscInt iva = 0; iva < ModelAData::NV; iva++) {
              V[iva] = rdtg * ModelARndm->variance1() / f;
              A[iva] = rdtg * ModelARndm->variance1() / f;
            }
            O4AlgebraHelper::O4Rotation(V, A, phi_n.f);
          } else {
            // This is the normal model A updatae
            for (PetscInt l = 0; l < ModelAData::Nphi; l++) {
              phi_o.f[l] = phi[k][j][i].f[l]; // Fill up the old values
              phi_n.f[l] = phi_o.f[l] + rdtg * ModelARndm->variance1();
            }
          }

          for (PetscInt l = 0; l < ModelAData::Nphi; l++) {
            heff.f[l] = (phi[k][j][i + 1].f[l] + phi[k][j][i - 1].f[l]) +
                        (phi[k][j + 1][i].f[l] + phi[k][j - 1][i].f[l]) +
                        (phi[k + 1][j][i].f[l] + phi[k - 1][j][i].f[l]);

            s_o += pow(phi_o.f[l], 2);
            s_n += pow(phi_n.f[l], 2);

            hphi_o += heff.f[l] * phi_o.f[l];
            hphi_n += heff.f[l] * phi_n.f[l];
          }
          hphi_o += (H * phi_o.f[0]); // Take care of the true H
          hphi_n += (H * phi_n.f[0]);

          double dS = -(hphi_n - hphi_o) + Lambda2 * (s_n - s_o) +
                      lambda / 4. * (s_n * s_n - s_o * s_o);

          // Downward step
          if (dS < 0) {
            for (PetscInt l = 0; l < ModelAData::Nphi; l++) {
              phinew[k][j][i].f[l] = phi_n.f[l];
            }
            monitor.increment_down(dS);
            continue;
          }
          // Process upward step
          double r = ModelARndm->uniform();
          if (r < exp(-dS)) {
            // keep the upward step w. probl exp(-dS)
            for (PetscInt l = 0; l < ModelAData::Nphi; l++) {
              phinew[k][j][i].f[l] = phi_n.f[l];
            }
            monitor.increment_up_yes(dS);
            continue;
          } else {
            // dphi is rejected. Set the action change dS=0
            monitor.increment_up_no(0.);
            continue;
          }
        }
      }
    }
    PetscLogEventEnd(loop, 0, 0, 0, 0);
    // Retstore the arrays
    PetscCall(DMDAVecRestoreArray(model->domain, phi_local, &phi));
    PetscCall(DMDAVecRestoreArray(model->domain, model->solution, &phinew));
  }

  return true;
}
void EulerLangevinHB::finalize() {
  VecDestroy(&phi_local);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    monitor.print(stdout);
  }
}

/////////////////////////////////////////////////////////////////////////

// The six possible faces of the cube on the boundary of cell A and cell B. For
// instance,  face_cases[0] (the first row), is an even site site for cell B,
// eoA=0, cell B is shifted by plus one in the x direction.  and not shited in
// the y and z directions. face_cases[5] (the last row) is odd (eoA=1) for cell
// A, and cell B is shifted in z by one unit
/* clang-format off */
g_face_case g_face_cases[3][2] = {
    0, 1, 0, 0,
    1, 1, 0, 0,
    0, 0, 1, 0,
    1, 0, 1, 0,
    0, 0, 0, 1,
    1, 0, 0, 1,
};
/* clang-format on */

// Update a pair of cells A and B with the Heat Bath
double modelg_update_charge_pair(const double &chi, const double &rms,
                                 const double &nA, const double &nB,
                                 o4_stepper_monitor &monitor) {
  double q = rms * ModelARndm->variance1();
  double dS =
      (pow(nA - q, 2) + pow(nB + q, 2) - pow(nA, 2) - pow(nB, 2)) / (2. * chi);
  // Downward step
  if (dS < 0) {
    monitor.increment_down(dS);
    return q;
  }
  // Process upward step
  double r = ModelARndm->uniform();
  if (r < exp(-dS)) {
    // keep the upward step w. probl exp(-dS)
    monitor.increment_up_yes(dS);
    return q;
  } else {
    // q is rejected. Set the action change dS=0
    monitor.increment_up_no(0.);
    return 0.;
  }
}

ModelGChargeHB::ModelGChargeHB(ModelA &in)
    : model(&in), qmonitor("Charge HB Steps") {
  DMCreateLocalVector(model->domain, &phi_local);
  DMCreateLocalVector(model->domain, &dn_local);
}

bool ModelGChargeHB::step(const double &dt) {

  // Get the ranges
  PetscInt ixs, iys, izs, nx, ny, nz;
  PetscCall(DMDAGetCorners(model->domain, &ixs, &iys, &izs, &nx, &ny, &nz));

  const auto &coeff = model->data.acoefficients;
  const PetscReal rdtsigma = sqrt(2. * dt * coeff.sigma());
  const PetscReal chi = coeff.chi;

  // Shuffle the order of xyz and the order of even and odd
  // to eliminate potential bias.
  std::array<int, 3> orderxyz{0, 1, 2};
  std::array<int, 2> ordereo{0, 1};
  if (model->rank == 0) {
    std::shuffle(orderxyz.begin(), orderxyz.end(), ModelARndm->generator());
    std::shuffle(ordereo.begin(), ordereo.end(), ModelARndm->generator());
  }
  MPI_Bcast(orderxyz.data(), 3, MPI_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(ordereo.data(), 2, MPI_INT, 0, PETSC_COMM_WORLD);

  // Checkerboard order ieo = even and odd sites
  for (int ixyz = 0; ixyz < 3; ixyz++) {
    for (int ieo = 0; ieo < 2; ieo++) {
      // get the face case that we will update
      g_face_case &face = g_face_cases[orderxyz[ixyz]][ordereo[ieo]];

      // take the global vector U and distribute to the local vector localU
      PetscCall(DMGlobalToLocalBegin(model->domain, model->solution,
                                     INSERT_VALUES, phi_local));
      PetscCall(DMGlobalToLocalEnd(model->domain, model->solution,
                                   INSERT_VALUES, phi_local));

      // Zero out the differences
      PetscCall(VecSet(dn_local, 0.));

      // Get the array to store the charge transfers
      data_node ***dn;
      PetscCall(DMDAVecGetArray(model->domain, dn_local, &dn));

      // Get pointer to local array
      data_node ***phi;
      PetscCall(DMDAVecGetArray(model->domain, phi_local, &phi));

      for (int k = izs; k < izs + nz; k++) {
        for (int j = iys; j < iys + ny; j++) {
          for (int i = ixs; i < ixs + nx; i++) {
            if ((k + j + i) % 2 != face.eoA) {
              continue;
            }
            for (int L = ModelAData::Nphi; L < ModelAData::Ndof; L++) {
              int iB = i + face.iB;
              int jB = j + face.jB;
              int kB = k + face.kB;
              const PetscScalar &nA = phi[k][j][i].x[L];
              const PetscScalar &nB = phi[kB][jB][iB].x[L];
              PetscScalar q =
                  modelg_update_charge_pair(chi, rdtsigma, nA, nB, qmonitor);

              dn[k][j][i].x[L] = -q;
              dn[kB][jB][iB].x[L] = q;
            }
          }
        }
      }
      PetscCall(DMDAVecRestoreArray(model->domain, dn_local, &dn));
      PetscCall(DMDAVecRestoreArray(model->domain, phi_local, &phi));
      PetscCall(DMLocalToGlobal(model->domain, dn_local, ADD_VALUES,
                                model->solution));
    }
  }

  return true;
};

void ModelGChargeHB::finalize() {
  PetscCallVoid(VecDestroy(&phi_local));
  PetscCallVoid(VecDestroy(&dn_local));

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    qmonitor.print(stdout);
  }
}

/////////////////////////////////////////////////////////////////////////

ModelGDiffusionStep::ModelGDiffusionStep(ModelA &in, const bool &implicit_step)
    : model(&in), use_implicit_step(implicit_step) {

  VecDuplicate(model->solution, &rhs);
  VecDuplicate(model->solution, &dn);
  DMCreateMatrix(model->domain, &J);

  double hx = model->data.hX();
  double hy = model->data.hY();
  double hz = model->data.hZ();
  double GammaOverD =
      model->data.acoefficients.Gamma() / model->data.acoefficients.D();
  Form3PointLaplacian(model->domain, J, hx, hy, hz, GammaOverD);

  MatConvert(J, MATSAME, MAT_INITIAL_MATRIX, &A);
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetFromOptions(ksp);
}

void ModelGDiffusionStep::finalize() {
  KSPDestroy(&ksp);
  MatDestroy(&A);
  MatDestroy(&J);
  VecDestroy(&dn);
  VecDestroy(&rhs);
}

// We are solving the diffusion equation for model G:
//
// \partial_t n = D nabla^2 n
//
// The fields phi,  in superfluid mode,  are normalized to f2,  so after
// updating the phi fields according the diffusion equation, we need to
// normalize the first four components of the solution to f2.
//
// If use_implicit method is true  the  second order Crank-Nicolson scheme is
// used. At first order Crank-Nicolson scheme we have
//
// (1/dtD  + J) n_+ = n/dtD
//
// Here J = - nabla^2,  dtD = dt * D and n_+ are the fields at time t + dt.
//
// At second order Crank-Nicolson scheme we have, similarly,
//
//  (2/dtD + J) n_+ = 2 n/dtD - J n
//
// If use_implicit_step is false,  we are using the explicit second order
// Runge-Kutta scheme.  The timestep should be smaller than 1/(12*D) to
// ensure stability.
bool ModelGDiffusionStep::step(const double &dt) {

  // Parameters needed
  const auto &coeff = model->data.acoefficients;
  const PetscReal &dtD = dt * coeff.D();
  const auto &f2 = model->data.f2();
  bool superfluidmode = model->data.ahandler.superfluidmode;

  if (use_implicit_step) {
    // This was used for the first order Crank-Nicolson scheme
    // VecCopy(model->solution, rhs);
    // VecScale(rhs, 1. / dtD);

    // Second order Crank-Nicolson scheme
    MatMult(J, model->solution, dn);
    VecCopy(model->solution, rhs);
    VecScale(rhs, 2. / dtD);
    VecAXPY(rhs, -1.0, dn); // rhs = rhs - J * n

    // This is for the second order Crank-Nicolson scheme
    MatCopy(J, A, SAME_NONZERO_PATTERN);
    // 2/dtD is for second order Crank-Nicolson scheme
    MatShift(A, 2. / dtD);

    // Actually solve the linear system
    KSPSetOperators(ksp, A, A);
    KSPSolve(ksp, rhs, model->solution);
  } else {
    // Explicit second order Runge-Kutta scheme
    MatMult(J, model->solution, dn);
    VecCopy(model->solution, rhs);
    VecAXPY(rhs, -0.5 * dtD, dn); // comput n at t + dt/2
    MatMult(J, rhs, dn);
    VecAXPY(model->solution, -1. * dtD, dn); //
  }

  if (superfluidmode) {
    // Normalize the first four components of the solution
    PetscInt i, j, k, xstart, ystart, zstart, xdimension, ydimension,
        zdimension;

    DMDAGetCorners(model->domain, &xstart, &ystart, &zstart, &xdimension,
                   &ydimension, &zdimension);

    G_node ***phi;
    PetscCall(DMDAVecGetArray(model->domain, model->solution, &phi));
    double f = sqrt(f2);
    for (k = zstart; k < zstart + zdimension; k++) {
      for (j = ystart; j < ystart + ydimension; j++) {
        for (i = xstart; i < xstart + xdimension; i++) {
          G_node::normalize_phi(phi[k][j][i].f, f);
        }
      }
    }
    PetscCall(DMDAVecRestoreArray(model->domain, model->solution, &phi));
  }

  return true;
}

// Form the Jacobian the Jabian generic interface
PetscErrorCode
ModelGDiffusionStep::Form3PointLaplacian(DM da, Mat J, const double &hx,
                                         const double &hy, const double &hz,
                                         const double &GammaOverD) {
  // Get the local information and store in info
  DMDALocalInfo info;
  DMDAGetLocalInfo(da, &info);
  PetscInt i, j, k, l;
  for (k = info.zs; k < info.zs + info.zm; k++) {
    for (j = info.ys; j < info.ys + info.ym; j++) {
      for (i = info.xs; i < info.xs + info.xm; i++) {
        for (l = 0; l < ModelAData::Ndof; l++) {
          PetscScalar r = (l < ModelAData::Nphi) ? GammaOverD : 1.0;
          // we define the column
          PetscInt nc = 0;
          MatStencil row, column[10];
          PetscScalar value[10];
          // here we insert the position of the row
          row.i = i;
          row.j = j;
          row.k = k;
          row.c = l;
          // here we define de position of the non-vansih column for the given
          // row in total there are 7*4 entries and nc is the total number of
          // column per row x direction
          column[nc].i = i - 1;
          column[nc].j = j;
          column[nc].k = k;
          column[nc].c = l;
          value[nc++] = -1. * r / (hx * hx);
          column[nc].i = i + 1;
          column[nc].j = j;
          column[nc].k = k;
          column[nc].c = l;
          value[nc++] = -1. * r / (hx * hx);
          // y direction
          column[nc].i = i;
          column[nc].j = j - 1;
          column[nc].k = k;
          column[nc].c = l;
          value[nc++] = -1. * r / (hy * hy);
          column[nc].i = i;
          column[nc].j = j + 1;
          column[nc].k = k;
          column[nc].c = l;
          value[nc++] = -1. * r / (hy * hy);
          // z direction
          column[nc].i = i;
          column[nc].j = j;
          column[nc].k = k - 1;
          column[nc].c = l;
          value[nc++] = -1. * r / (hz * hz);
          column[nc].i = i;
          column[nc].j = j;
          column[nc].k = k + 1;
          column[nc].c = l;
          value[nc++] = -1. * r / (hz * hz);

          // The central element need a loop over the flavour index of the
          // column (is a full matrix in the flavour index )
          column[nc].i = i;
          column[nc].j = j;
          column[nc].k = k;
          column[nc].c = l;
          value[nc++] =
              2.0 * r / (hx * hx) + 2.0 * r / (hy * hy) + 2.0 * r / (hz * hz);

          // here we set the matrix. Petsc wraps the boundary conditions
          // autmatically
          MatSetValuesStencil(J, 1, &row, nc, column, value, INSERT_VALUES);
        }
      }
    }
  }
  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  return (0);
}

/////////////////////////////////////////////////////////////////////////

bool ModelGExplicitDiffusionStep::evolveLocalSolution(const double &dt,
                                                      G_node ***phi,
                                                      G_node ***phinew){

  auto &data = model->data;

  // Get a local vector with ghost cells
  DM da = model->domain;

  const auto &coeff = data.acoefficients;
  PetscReal H[4] = {coeff.H, 0., 0., 0.};
  if (data.ahandler.superfluidmode) {
    H[0] *= coeff.sigmabyf(data.atime.t());
  }
  const PetscReal axx = pow(1. / data.hX(), 2);
  const PetscReal ayy = pow(1. / data.hY(), 2);
  const PetscReal azz = pow(1. / data.hZ(), 2);

  PetscReal dtG = dt * coeff.Gamma();
  PetscReal dtD = dt * coeff.D();
  PetscReal f = sqrt(data.f2());

  PetscInt xstart, ystart, zstart, xdimension, ydimension, zdimension;
  DMDAGetCorners(da, &xstart, &ystart, &zstart, &xdimension, &ydimension,
                 &zdimension);

  // Loop over central elements
  for (PetscInt k = zstart; k < zstart + zdimension; k++) {
    for (PetscInt j = ystart; j < ystart + ydimension; j++) {
      for (PetscInt i = xstart; i < xstart + xdimension; i++) {
        // First evolve the momenta nab

        G_node &centralPhi = phi[k][j][i];
        G_node &phixplus = phi[k][j][i + 1];
        G_node &phixminus = phi[k][j][i - 1];
        G_node &phiyplus = phi[k][j + 1][i];
        G_node &phiyminus = phi[k][j - 1][i];
        G_node &phizplus = phi[k + 1][j][i];
        G_node &phizminus = phi[k - 1][j][i];

        for (PetscInt l = 0; l < ModelAData::Nphi; l++) {

          // Compute the nabla^2 phi
          phinew[k][j][i].f[l] +=
              dtG * (axx * (phixplus.f[l] + phixminus.f[l]) +
                     ayy * (phiyplus.f[l] + phiyminus.f[l]) +
                     azz * (phizplus.f[l] + phizminus.f[l]) -
                     2. * (axx + ayy + azz) * centralPhi.f[l]);

          // Add the mass term (H[l]  is  H * sigma/f)
          phinew[k][j][i].f[l] += dtG * H[l];
        }
        if (data.ahandler.superfluidmode) {
          // In superfluid mode we need to normalize the first four components
          // of the solution to f2
          G_node::normalize_phi(phinew[k][j][i].f, f);
        }

        for (PetscInt l = 0; l < ModelAData::NA; l++) {
          // Compute the nabla^2 A
          phinew[k][j][i].A[l] +=
              dtD * (axx * (phixplus.A[l] + phixminus.A[l]) +
                     ayy * (phiyplus.A[l] + phiyminus.A[l]) +
                     azz * (phizplus.A[l] + phizminus.A[l]) -
                     2. * (axx + ayy + azz) * centralPhi.A[l]);
        }

        for (PetscInt l = 0; l < ModelAData::NV; l++) {
          phinew[k][j][i].V[l] +=
              dtD * (axx * (phixplus.V[l] + phixminus.V[l]) +
                     ayy * (phiyplus.V[l] + phiyminus.V[l]) +
                     azz * (phizplus.V[l] + phizminus.V[l]) -
                     2. * (axx + ayy + azz) * centralPhi.V[l]);
        }
      }
    }
  }

  return true;
}

// call to do a diffusion step that evolves current solution and updates it 
bool ModelGExplicitDiffusionStep::step(const double &dt) {

  auto &data = model->data;

  // Get a local vector (buffer) with ghost cells
  DM da = model->domain;
  Vec localU;
  PetscCall(DMGetLocalVector(da, &localU));

  // Fill in the ghost celss with mpicalls
  PetscCall(DMGlobalToLocalBegin(da, model->solution, INSERT_VALUES, localU));
  PetscCall(DMGlobalToLocalEnd(da, model->solution, INSERT_VALUES, localU));

  G_node ***phi, ***phinew;
  // get access to arrays indexed using the local dimensions
  PetscCall(DMDAVecGetArrayRead(da, localU, &phi));
  // this one needs to also have write access & points to the solution directly
  // but careful! this only works since we only write to phinew and only read from it locally
  // otherwise, corruption through ghost cells is possible
  PetscCall(DMDAVecGetArray(da, model->solution, &phinew));

  // call subroutine that computes phinew
  evolveLocalSolution(dt, phi, phinew);

  // restore pointer arrays
  PetscCall(DMDAVecRestoreArray(da, model->solution, &phinew));
  PetscCall(DMDAVecRestoreArrayRead(da, localU, &phi));

  // restore local vectors
  PetscCall(DMRestoreLocalVector(da, &localU));

  return true;
}

// call to do a diffusion step that evolves current solution which was deep copied into
// solution_coarsened nsteps coarsen-steps (each = 3 diffusion steps) but does not
// update the model solution, but stores it into input solution_coarsened
bool ModelGExplicitDiffusionStep::step_coarsening(const double &dt,
                                                  Vec *solution_coarsened,
                                                  int nsteps) {

  // each coarsen step consists of 3 diffusion steps with dt/3
  int ndiffusion_steps = 3 * nsteps;

  // Get a local vector with ghost cells
  DM da = model->domain;
  Vec localU;
  PetscCall(DMGetLocalVector(da, &localU));

  // pointer arrays
  G_node ***phi, ***phinew;

  for (int i = 0; i < ndiffusion_steps; i++) {
    // Fill in the ghost cells with mpicalls
    PetscCall(DMGlobalToLocalBegin(da, *solution_coarsened, INSERT_VALUES, localU));
    PetscCall(DMGlobalToLocalEnd(da, *solution_coarsened, INSERT_VALUES, localU));
    // all mpi calls complete: localU is up to date

    // get access to arrays indexed using the local dimensions
    PetscCall(DMDAVecGetArrayRead(da, localU, &phi));
    // this one needs to also have write access & points to solution_coarsened directly
    // but careful! this only works since we only write to phinew and only read from it locally
    // otherwise, corruption through ghost cells is possible
    PetscCall(DMDAVecGetArray(da, *solution_coarsened, &phinew));
  
    // call subroutine that computes phinew
    evolveLocalSolution(dt/3., phi, phinew);

    // restore localUnew from phinew
    PetscCall(DMDAVecRestoreArray(da, *solution_coarsened, &phinew));
    PetscCall(DMDAVecRestoreArrayRead(da, localU, &phi));
  }

  // Restore local vectors after loop
  PetscCall(DMRestoreLocalVector(da, &localU));

  return true;
}

/////////////////////////////////////////////////////////////////////////

PV2HBSplit::PV2HBSplit(ModelA &in, const std::string &insteps,
                       const bool &ideal, const bool &heatbath,
                       const bool &diffusion)
    : model(&in), pv2(in), hbPhi(in), hbN(in), steps(insteps),
      include_ideal(ideal), include_heatbath(heatbath),
      include_diffusion(diffusion), A_count(0), B_count(0), C_count(0) {

  // Check if string contains only A, B, C characters
  for (char c : steps) {
    if (c != 'A' and c != 'B' and c != 'C') {
      throw std::invalid_argument(
          "String contains invalid character: " + std::string(1, c) +
          ". Only A, B, C are allowed.");
    }
  }
  // Count occurrences of each character
  for (char c : steps) {
    switch (c) {
    case 'A':
      ++A_count;
      break;
    case 'B':
      ++B_count;
      break;
    case 'C':
      ++C_count;
      break;
    }
  }
}

bool PV2HBSplit::step(const double &dt) {

  // Compute the ideal step with the position verlet integrator.
  // Record the energy for the accept reject step

  PetscLogEvent ideal_log, hb_log, qhb_log;

  // Format of steps is ABBBABBBC
  PetscLogEventRegister("IdealStep", 0, &ideal_log);
  PetscLogEventRegister("HBStep", 0, &hb_log);
  PetscLogEventRegister("QHBStep", 0, &qhb_log);

  for (char s : steps) {
    if (s == 'A' and include_ideal) {
      PetscLogEventBegin(ideal_log, 0, 0, 0, 0);
      pv2.step(dt / A_count);
      PetscLogEventEnd(ideal_log, 0, 0, 0, 0);
    } else if (s == 'B' and include_heatbath) {
      PetscLogEventBegin(hb_log, 0, 0, 0, 0);
      hbPhi.step(dt / B_count);
    } else if (s == 'C' and include_diffusion) {
      PetscLogEventBegin(qhb_log, 0, 0, 0, 0);
      hbN.step(dt / C_count);
      PetscLogEventEnd(qhb_log, 0, 0, 0, 0);
    }
  }
  return true;
}
