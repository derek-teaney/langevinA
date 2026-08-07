#ifndef MEASURER
#define MEASURER

#include "ModelA.h"
// #include "make_unique.h"
#include "nvector.h"
#include "Stepper.h"
#include <array>
#include <complex>
#include <fftw3.h>
#include <vector>

////////////////////////////////////////////////////////////////////////
// Class to compute the fourier transform of the 1d array. The fourier transform
// is defined as 1/N \sum_x e^{ikx} W(x) = W(k) where W(x) is the array of size
// N and W(k) is the complex array of size N/2 + 1
class measurer_fft {

private:
  // Create the storage space for the plan
  std::vector<double> in_ptr;
  std::vector<std::complex<double>> out_ptr;
  fftw_plan plan;

public:
  measurer_fft(const size_t N) : in_ptr(N, 0.), out_ptr(N / 2 + 1, 0.) {
    // Create the plane for the fourier transform
    plan = fftw_plan_dft_r2c_1d(
        N, in_ptr.data(), reinterpret_cast<fftw_complex *>(out_ptr.data()),
        FFTW_MEASURE);
  }
  ~measurer_fft() { fftw_destroy_plan(plan); }

  void execute_row(double *in, std::complex<double> *out) {
    for (size_t i = 0; i < in_ptr.size(); i++) {
      in_ptr[i] = in[i];
    }
    fftw_execute(plan);
    for (size_t i = 0; i < out_ptr.size(); i++) {
      out[i] = out_ptr[i] / in_ptr.size();
    }
  }

  // Compute the fourier transform of an array in
  void execute(nvector<double, 2> &in, nvector<std::complex<double>, 2> &out) {
    for (size_t i = 0; i < in.N[0]; i++) {
      execute_row(&in(i, 0), &out(i, 0));
    }
  }
};

////////////////////////////////////////////////////////////////////////
// Class to compute a real-to-complex 3D fourier transform of a scalar field
// on a periodic box.  The transform uses the same FFTW backend and the same
// normalization convention as measurer_fft above, generalized to 3D:
//
//   W(k) = (1 / (n0*n1*n2)) \sum_x e^{-i k.x} W(x)
//
// The dimensions are (n0, n1, n2) in row-major (C) order, so that the last
// index n2 is contiguous in memory.  The real input has size n0*n1*n2 and the
// complex output (FFTW r2c half spectrum) has size n0*n1*(n2/2+1).
class measurer_fft_3d {

private:
  size_t n0, n1, n2;
  std::vector<double> in_ptr;
  std::vector<std::complex<double>> out_ptr;
  fftw_plan plan;

public:
  measurer_fft_3d(const size_t N0, const size_t N1, const size_t N2)
      : n0(N0), n1(N1), n2(N2), in_ptr(N0 * N1 * N2, 0.),
        out_ptr(N0 * N1 * (N2 / 2 + 1), 0.) {
    plan = fftw_plan_dft_r2c_3d(
        N0, N1, N2, in_ptr.data(),
        reinterpret_cast<fftw_complex *>(out_ptr.data()), FFTW_MEASURE);
  }
  ~measurer_fft_3d() { fftw_destroy_plan(plan); }

  // Transform the real volume in (size n0*n1*n2) into the complex half
  // spectrum out (size n0*n1*(n2/2+1)), applying the 1/(n0*n1*n2) normalization.
  void execute(const std::vector<double> &in,
               std::vector<std::complex<double>> &out) {
    std::copy(in.begin(), in.end(), in_ptr.begin());
    fftw_execute(plan);
    const double norm = static_cast<double>(n0 * n1 * n2);
    for (size_t i = 0; i < out_ptr.size(); i++) {
      out[i] = out_ptr[i] / norm;
    }
  }
};

////////////////////////////////////////////////////////////////////////
// Computes slices of the fields and their fourier transforms
//
// The fields are labelled by U = phi, A, V, and phi2 and hence has dimensions
// NObs = Nphi +  NA +  NV  + 1 respectively.  The field phi2 is the square of
// the  field phi_a.
//
// Our fourier transform conventions are  (1/N_x) \sum_x e^{ikx} W(x) = W(k)
// where W(x)  is the average of the fields over the y and z directions, i.e.
// W(x) = 1/N_y * 1/N_z \sum_{y,z} U(x,y,z)
//
// The array wallX has dimensions NObs x N and contains the slices of the fields
// in the x direction. The array wallX_k has dimensions NObs x N/2 + 1 and
// contains the fourier transform of the slices of the fields in the x
// direction.
//
// The fields are rotated in the direction of the vev and the rotated fields are
// stored in the array wallX_rotated. The array wallXPhase has dimensions
// NObsPhase x N and contains the slices of the fields in the x direction of the
// rotated fields.  Finally the array wallXPhase_k has dimensions NObsPhase x
// N/2 + 1 and contains the fourier transform of the slices of the fields in the
// x direction of the rotated fields.
//
// The rotated fields are define as follows: The zero mode of the field phi
// defines a unit four vector n_a. The rotated fields are defined as follows:
//
// sigma = n_a phi_a,
// pi_b = phi_a - sigma n_b,
// A_a = rho_{ab} n_b,
// V_a = rhotilde_{ab} n_b,
// phi2.
//
// Here rho_{ab} is the field A_a and V_a are the
// axial and vector densities respectively and rhotilde_{ab} is the dual of the
// field rho_{ab}, i.e. multiplied by the epsilon tensor epsilon_{abcd} which
// swaps the vector and axial vector pieces.  Thus the dimension of the rotated
// fields is NObsPhase = 3*Nphi + 2 where Nphi is the number of scalar fields.
class Measurer {
public:
  // Arrays of size Nobs contain X=(phi[1..Nphi], q[1...Nq], phi2)
  static const PetscInt NObs =
      ModelAData::Nphi + ModelAData::NA + ModelAData::NV + 1;

  // Scalar data this is k=0 mode by itself
  static const PetscInt NScalars = NObs;
  std::vector<PetscScalar> OAverage;

  // energy data (different parts of H)
  static const PetscInt NEnergy = 5;
  std::vector<PetscScalar> Energy;
  std::vector<PetscScalar> EnergyRotated;
  std::vector<PetscScalar> EnergyPhase;

  // Spherically-averaged Fourier readout of the topological charge density
  // q(x) of the O(4) field (see computeTopChargeFourier in measurer.cxx).
  //
  // NOTE: the spherical (azimuthal) average assumes an isotropic medium.  It
  // retains length-scale / ordering information but discards angular
  // (lattice-symmetry / orientation) information.
  int NtopchargeBins = 0;
  double topcharge_dk = 0.;             // radial bin width = 2*pi/L
  std::vector<double> topcharge_Sk;     // S(k) = <|qtilde|^2> per radial shell
  std::vector<double> topcharge_kbins;  // |k| bin centers, = m*dk
  std::vector<double> topcharge_Nshell; // number of modes per shell
  std::complex<double> topcharge_zero{0., 0.}; // qtilde(k=0), DC component

  // First dimension is NObs, last is spatial index x=0...N
  nvector<PetscScalar, 2> wallX;
  nvector<PetscScalar, 2> wallY;
  nvector<PetscScalar, 2> wallZ;

  // First dimension is NObs, last dimension is fourier index k=0..N/2+1
  nvector<std::complex<double>, 2> wallX_k;
  nvector<std::complex<double>, 2> wallY_k;
  nvector<std::complex<double>, 2> wallZ_k;

  // First dimension is NObsRotated, last dimension is fourier index.
  // Array of size NobsRotate contains XPhase = (sigma, pi[1..Nphi],
  // q[1..2*Nphi], phi2)
  static const PetscInt NObsRotated = 3 * ModelAData::Nphi + 2;
  nvector<std::complex<double>, 2> wallX_k_rotated;
  nvector<std::complex<double>, 2> wallY_k_rotated;
  nvector<std::complex<double>, 2> wallZ_k_rotated;

  // Array of size NobsPhase contains XPhase = (sigma, pi[1..Nphi],
  // q[1..2*Nphi], phi2)
  // The first dimension is NObsPhase, the second domenions is the spatial index
  static const PetscInt NObsPhase = 3 * ModelAData::Nphi + 2;
  nvector<PetscScalar, 2> wallXPhase;
  nvector<PetscScalar, 2> wallYPhase;
  nvector<PetscScalar, 2> wallZPhase;

  // First dimension is NObsPhase, last dimension is the fourier index 0..N/2+1
  nvector<std::complex<double>, 2> wallXPhase_k;
  nvector<std::complex<double>, 2> wallYPhase_k;
  nvector<std::complex<double>, 2> wallZPhase_k;

  // Array of size NObsCoarse contains XCoarse = (sigma, pi[1..Nphi],
  // q[1..2*Nphi], phi2)
  // The first dimension is NObsCoarse, the second dimension is the spatial index
  static const PetscInt NObsCoarse = 3 * ModelAData::Nphi + 2;
  Vec solution_coarsened;
  // Selected coarsen levels to record: ncoarsen_start, ncoarsen_start+stride, ...
  std::vector<int> coarsen_levels;
  // Vectors of length coarsen_levels.size(); index c holds result for coarsen_levels[c] steps
  std::vector<nvector<PetscScalar, 2>> wallXCoarse;
  std::vector<nvector<PetscScalar, 2>> wallYCoarse;
  std::vector<nvector<PetscScalar, 2>> wallZCoarse;

  // First dimension is NObsCoarse, last dimension is the fourier index 0..N/2+1
  std::vector<nvector<std::complex<double>, 2>> wallXCoarse_k;
  std::vector<nvector<std::complex<double>, 2>> wallYCoarse_k;
  std::vector<nvector<std::complex<double>, 2>> wallZCoarse_k;

public:
  Measurer(ModelA *ptr) : model(ptr) {
    N = model->data.NX;
    if (model->data.NX != model->data.NY || model->data.NX != model->data.NZ) {
      PetscPrintf(
          PETSC_COMM_WORLD,
          "Nx, Ny, and Nz must be equal for the correlation analysis to work");
      throw(
          "Nx, Ny, and Nz must be equal for the correlation analysis to work");
    }

    // Only need the fft for derived observables
    int rank = -1;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) {
       fftw = make_unique<measurer_fft>(N);
    }


    wallX.resize(NObs, N);
    wallY.resize(NObs, N);
    wallZ.resize(NObs, N);

    wallX_k.resize(NObs, N / 2 + 1);
    wallY_k.resize(NObs, N / 2 + 1);
    wallZ_k.resize(NObs, N / 2 + 1);

    wallX_k_rotated.resize(NObsRotated, N / 2 + 1);
    wallY_k_rotated.resize(NObsRotated, N / 2 + 1);
    wallZ_k_rotated.resize(NObsRotated, N / 2 + 1);

    wallXPhase.resize(NObsPhase, N);
    wallYPhase.resize(NObsPhase, N);
    wallZPhase.resize(NObsPhase, N);

    wallXPhase_k.resize(NObsPhase, N / 2 + 1);
    wallYPhase_k.resize(NObsPhase, N / 2 + 1);
    wallZPhase_k.resize(NObsPhase, N / 2 + 1);

    {
      const auto &h = model->data.ahandler;
      for (int s = h.ncoarsen_start; s <= h.ncoarsen_steps; s += h.ncoarsen_stride)
        coarsen_levels.push_back(s);
    }
    int nout = static_cast<int>(coarsen_levels.size());
    wallXCoarse.resize(nout);
    wallYCoarse.resize(nout);
    wallZCoarse.resize(nout);
    wallXCoarse_k.resize(nout);
    wallYCoarse_k.resize(nout);
    wallZCoarse_k.resize(nout);
    for (int c = 0; c < nout; c++) {
      wallXCoarse[c].resize(NObsCoarse, N);
      wallYCoarse[c].resize(NObsCoarse, N);
      wallZCoarse[c].resize(NObsCoarse, N);
      wallXCoarse_k[c].resize(NObsCoarse, N / 2 + 1);
      wallYCoarse_k[c].resize(NObsCoarse, N / 2 + 1);
      wallZCoarse_k[c].resize(NObsCoarse, N / 2 + 1);
    }

    // Set up the radial binning for the spherically averaged topological
    // charge spectrum.  The bin width is the reciprocal-lattice spacing
    // dk = 2*pi/L (LX == LY == LZ here since NX == NY == NZ).  The largest
    // possible |k| integer magnitude is sqrt(3)*(N/2), so we need
    // floor(sqrt(3)*N/2 + 0.5)+1 bins with center |k|_m = m*dk.
    topcharge_dk = 2.0 * M_PI / model->data.LX;
    const double rmax = sqrt(3.0) * (0.5 * static_cast<double>(N));
    NtopchargeBins = static_cast<int>(floor(rmax + 0.5)) + 1;
    topcharge_Sk.assign(NtopchargeBins, 0.);
    topcharge_kbins.assign(NtopchargeBins, 0.);
    topcharge_Nshell.assign(NtopchargeBins, 0.);
    for (int m = 0; m < NtopchargeBins; m++) {
      topcharge_kbins[m] = m * topcharge_dk;
    }
    // Only rank 0 assembles the full q(x) volume and takes its 3D FFT.
    if (rank == 0) {
      fftw3d = make_unique<measurer_fft_3d>(N, N, N);
    }

    // create unique diffusion stepper
    diffuser = std::make_unique<ModelGExplicitDiffusionStep>(*model);
    // create global vector that stores coarsened solution
    DMCreateGlobalVector(model->domain, &solution_coarsened);

    PetscLogEventRegister("TopCharge_measure", 0, &topcharge_log);
    PetscLogEventRegister("Energy_measure", 0, &energy_log);
    PetscLogEventRegister("Normal+Phase_measure", 0, &convt_log);
    PetscLogEventRegister("Coarsen_measure", 0, &coarsen_log);
    PetscLogEventRegister("Derived_measure", 0, &derived_log);
  }

  virtual ~Measurer() {}

  // Takes the vector, solution,  and computes the walls and their fourier
  // transforms and the rotated versions.  On output the arrays wallX, wallX_k,
  // wallX_rotated, as well as Y and Z are filled with the data. This can be
  // accessed to write the data to disk.
  void measure(Vec *solution) {

    int rank = -1;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscLogEventBegin(convt_log, 0, 0, 0, 0);
    computeSliceAverage(solution);
    computeSliceAveragePhase(solution);
    PetscLogEventEnd(convt_log, 0, 0, 0, 0);
    
    PetscLogEventBegin(energy_log, 0, 0, 0, 0);
    computeEnergy();
    computeEnergyRotated();
    computeEnergyPhase();
    PetscLogEventEnd(energy_log, 0, 0, 0, 0);

    // Take the FFT and other steps based on the data collected
    if (rank == 0) {
      PetscLogEventBegin(derived_log, 0, 0, 0, 0);
      computeDerivedObs();
      PetscLogEventEnd(derived_log, 0, 0, 0, 0);
    }
  }

  void measure_coarsen(Vec * solution){
    int rank = -1;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    PetscLogEventBegin(coarsen_log, 0, 0, 0, 0);
    computeSliceAverageCoarsened(solution);
    PetscLogEventEnd(coarsen_log, 0, 0, 0, 0);
    
    // Take the FFT and other steps based on the data collected
    if (rank == 0) {
      PetscLogEventBegin(derived_log, 0, 0, 0, 0);
      computeDerivedObs_coarsen();
      PetscLogEventEnd(derived_log, 0, 0, 0, 0);
    }
  }

  // Computes the topological charge density q(x), its 3D FFT and the
  // spherically averaged power spectrum S(k) together with the zero mode.
  // computeTopChargeFourier is collective (it reduces q(x) to rank 0) so this
  // must be called on all ranks; the FFT/binning happens internally on rank 0.
  void measure_topcharge(Vec *solution) {
    PetscLogEventBegin(topcharge_log, 0, 0, 0, 0);
    computeTopChargeFourier(solution);
    PetscLogEventEnd(topcharge_log, 0, 0, 0, 0);
  }

  ModelA *getModel() { return model; }
  PetscInt getN() { return N; }
  int getNTopchargeBins() { return NtopchargeBins; }
  double getTopchargeDk() { return topcharge_dk; }
  int getNCoarsenOutputs() { return static_cast<int>(coarsen_levels.size()); }
  const std::vector<int> &getCoarsenLevels() { return coarsen_levels; }

private:
  void computeSliceAverage(Vec *solution);
  void computeSliceAveragePhase(Vec *solution);
  void computeSliceAverageCoarsened(Vec *solution);
  void computeEnergy();
  void computeEnergyRotated();
  void computeEnergyPhase();
  void computeDerivedObs();
  void computeDerivedObs_coarsen();
  void computeTopChargeFourier(Vec *solution);

  ModelA *model;
  PetscInt N;

  // FFT engine using the fftw3 library
  std::unique_ptr<measurer_fft> fftw;
  // 3D FFT engine for the topological charge density (rank 0 only)
  std::unique_ptr<measurer_fft_3d> fftw3d;
  // Diffusion stepper to coarsen/diffuse the solution
  std::unique_ptr<ModelGExplicitDiffusionStep> diffuser;

  PetscLogEvent energy_log, convt_log, coarsen_log, derived_log, topcharge_log;
};

#endif
