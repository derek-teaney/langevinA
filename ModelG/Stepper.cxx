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
    MPI_Bcast(&reject, 1, MPIU_BOOL, 0, MPI_COMM_WORLD);

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
  const PetscReal H[4] = {coeff.H, 0., 0., 0.};
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
        for (int s = 0; s < ModelAData::Nphi; ++s) {

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

  const auto &coeff = model->data.acoefficients;
  const double &Lambda2 = 0.5 * (2. * 3. + model->data.mass()); // mass term
  const double &lambda = coeff.lambda;
  const double &H = coeff.H;
  const PetscReal rdtg = sqrt(2. * dt * coeff.gamma);

  PetscLogEvent communication, random, loop;
  PetscLogEventRegister("Communication", 0, &communication);
  PetscLogEventRegister("Random", 0, &random);
  PetscLogEventRegister("Loop", 0, &loop);

  // Checkerboard order ieo = even and odd sites
  for (int ieo = 0; ieo < 2; ieo++) {
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
    for (int k = izs; k < izs + nz; k++) {
      for (int j = iys; j < iys + ny; j++) {
        for (int i = ixs; i < ixs + nx; i++) {
          if ((k + j + i) % 2 != ieo) {
            continue;
          }

          PetscScalar hphi_o = 0.; // old and new h*phi
          PetscScalar hphi_n = 0.; // old and new h*phi
          PetscScalar s_o = 0.;    // old sum of phi**2
          PetscScalar s_n = 0.;    // new sum of phi**2

          for (int l = 0; l < ModelAData::Nphi; l++) {
            heff.f[l] = (phi[k][j][i + 1].f[l] + phi[k][j][i - 1].f[l]) +
                        (phi[k][j + 1][i].f[l] + phi[k][j - 1][i].f[l]) +
                        (phi[k + 1][j][i].f[l] + phi[k - 1][j][i].f[l]);

            phi_o.f[l] = phi[k][j][i].f[l];
            phi_n.f[l] = phi_o.f[l] + rdtg * ModelARndm->variance1();

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
            for (int l = 0; l < ModelAData::Nphi; l++) {
              phinew[k][j][i].f[l] = phi_n.f[l];
            }
            monitor.increment_down(dS);
            continue;
          }
          // Process upward step
          double r = ModelARndm->uniform();
          if (r < exp(-dS)) {
            // keep the upward step w. probl exp(-dS)
            for (int l = 0; l < ModelAData::Nphi; l++) {
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

PV2HBSplit::PV2HBSplit(ModelA &in, const std::array<unsigned int, 2> &scounts,
                       const bool &nodiffuse, const bool &onlydiffuse)
    : model(&in), pv2(in), hbPhi(in), hbN(in), stepcounts(scounts),
      nodiffusion(nodiffuse), onlydiffusion(onlydiffuse) {}

bool PV2HBSplit::step(const double &dt) {

  // Compute the ideal step with the position verlet integrator.
  // Record the energy for the accept reject step

  PetscLogEvent ideal_log, hb_log, qhb_log;

  // Format of steps is ABBB,ABBB,C  for (3,2)
  if (onlydiffusion) {
    // pass
  } else {
    PetscLogEventRegister("IdealStep", 0, &ideal_log);
    PetscLogEventRegister("HBStep", 0, &hb_log);

    for (size_t i1 = 0; i1 < stepcounts[1]; i1++) {
      PetscLogEventBegin(ideal_log, 0, 0, 0, 0);
      pv2.step(dt / (stepcounts[1]));
      PetscLogEventEnd(ideal_log, 0, 0, 0, 0);

      PetscLogEventBegin(hb_log, 0, 0, 0, 0);
      for (size_t i0 = 0; i0 < stepcounts[0]; i0++) {
        hbPhi.step(dt / (stepcounts[0] * stepcounts[1]));
      }
      PetscLogEventEnd(hb_log, 0, 0, 0, 0);
    }
  }

  if (nodiffusion) {
    // pass
  } else {
    PetscLogEventRegister("QHBStep", 0, &qhb_log);

    PetscLogEventBegin(qhb_log, 0, 0, 0, 0);
    hbN.step(dt);
    PetscLogEventEnd(qhb_log, 0, 0, 0, 0);
  }

  return true;
}
