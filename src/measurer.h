#ifndef MEASURER
#define MEASURER

#include "ModelA.h"
// #include "make_unique.h"
#include "nvector.h"
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
  }

  virtual ~Measurer() {}

  // Takes the vector, solution,  and computes the walls and their fourier
  // transforms and the rotated versions.  On output the arrays wallX, wallX_k,
  // wallX_rotated, as well as Y and Z are filled with the data. This can be
  // accessed to write the data to disk.
  void measure(Vec *solution) {
    int rank = -1;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    computeSliceAverage(solution);
    computeSliceAveragePhase(solution);

    // Take the FFT and other steps based on the data collected
    if (rank == 0) {
      computeDerivedObs();
    }
  }

  ModelA *getModel() { return model; }
  PetscInt getN() { return N; }

private:
  void computeSliceAverage(Vec *solution);
  void computeSliceAveragePhase(Vec *solution);
  void computeDerivedObs();

  ModelA *model;
  PetscInt N;

  // FFT engine using the fftw3 library
  std::unique_ptr<measurer_fft> fftw;
};

#endif
