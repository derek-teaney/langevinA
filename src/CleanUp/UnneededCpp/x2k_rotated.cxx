
#include "hdf5.h"
#include "ntuple.h"
#include <array>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <string>
#include <vector>

inline bool exists_test(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

inline std::array<std::complex<double>, 4>
contract_rho(const std::array<std::complex<double>, 4> &n,
             const std::array<std::complex<double>, 6> &rho) {
  return {-(n[1]*rho[0]) - n[2]*rho[1] - n[3]*rho[2],n[0]*rho[0] + n[3]*rho[4] - n[2]*rho[5], n[0]*rho[1] - n[3]*rho[3] + n[1]*rho[5],n[0]*rho[2] + n[2]*rho[3] - n[1]*rho[4]} ;
}

inline std::array<std::complex<double>, 6>
dualize_rho(const std::array<std::complex<double>, 6> &rho) {
  return std::array<std::complex<double>, 6>{rho[3], rho[4], rho[5],
                                             rho[0], rho[1], rho[2]};
}

// The purpose of the program is to take the fourier transform of a "wall"
// like dataset, and rotate it into the rotated frame.  Example:
//
// ./x2k.exe foo.h5 wally_k
//
// This will open foo.h5, and create foo_out.h5 if does not exist.  foo_out.h5
// will hall wally_k_rotated, which is the rotated version of the existing dataset.
// 
int main(int argc, char **argv) {

  if (argc < 3) {
    std::cout << "Usage x2k_rotated filename_out.h5 datasetname" << std::endl;
    exit(0);
  }

  std::string fin(argv[1]);
  std::string dsetin(argv[2]);

  std::string fout = fin;
  if (fout.rfind("_out.h5")==std::string::npos) {
    std::cout << "Usage x2k_rotated filename_out.h5 datasetname" << std::endl;
    exit(0);
  }

  std::string dsetout = dsetin + "_rotated";
  std::cout << "Creating dataset " << dsetout << " in " << fout << std::endl;


  // open the output file
  hid_t fileout_id;
  if (exists_test(fout)) {
    fileout_id = H5Fopen(fout.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  } else {
    std::cout << "Unable to open " << fout << std::endl;
  }
  
  // This is the input dset
  if (!H5Lexists(fileout_id, dsetin.c_str(), H5P_DEFAULT)) {
    std::cout << "The input dataset " << dsetin << " does not exist" << std::endl;
    H5Fclose(fileout_id);
  }

  ntuple<3> in(dsetin, fileout_id);
  // get the dimensions of the input tuple
  auto N = in.getN();

  // Remove the old output if it exists
  if (H5Lexists(fileout_id, dsetout.c_str(), H5P_DEFAULT)) {
    H5Ldelete(fileout_id, dsetout.c_str(), H5P_DEFAULT);
  }

  // Construct the dimensions of the output tuple.
  // 13 fields, lattice size, 0 and 1 for complex numbers
  std::array<size_t, 3> NK{13, N[1], N[2]};

  ntuple<3> out(NK, dsetout, fileout_id);

  std::array<std::complex<double>, 4> n{};
  std::array<std::complex<double>, 4> phi{};
  std::array<std::complex<double>, 4> phir{};
  std::array<std::complex<double>, 6> rho{};
  std::array<std::complex<double>, 4> V{};
  std::array<std::complex<double>, 4> A{};

  for (size_t i = 0; i < in.nrows(); i++) {
    in.readrow(i);

    double norm = 0.;
    for (size_t a = 0; a < 4; a++) {
      n[a] = {in.Row({a, 0, 0}), 0.};
      norm += std::norm(n[a]);
    }
    norm = sqrt(norm);
    if (norm < std::numeric_limits<double>::min()) {
      n = std::array<std::complex<double>, 4>{1, 0, 0, 0};
    } else {
      for (size_t a = 0; a < 4; a++) {
        n[a] /= norm;
      }
    }
    
    //testing
    //n = std::array<complex<double>, 4>{1, 0, 0, 0};

    for (size_t k = 0; k < N[1]; k++) {
      for (size_t a = 0; a < 4; a++) {
        phi[a] = {in.Row({a, k, 0}), in.Row({a, k, 1})};
      }

      for (size_t ab = 0; ab < 6; ab++) {
        rho[ab] = {in.Row({4 + ab, k, 0}),
                   in.Row({4 + ab, k, 1})};
      }

      std::complex<double> phis(0.);
      for (size_t a = 0; a < 4; a++) {
        phis += n[a] * phi[a];
      }
      for (size_t a = 0; a < 4; a++) {
        phir[a] = phi[a] - n[a] * phis;
      }

      A = contract_rho(n, rho);
      V = contract_rho(n, dualize_rho(rho));
      out.Row({0, k, 0}) = phis.real();
      out.Row({0, k, 1}) = phis.imag();
      for (size_t a = 0; a < 4; a++) {
        out.Row({1 + a, k, 0}) = phir[a].real();
        out.Row({1 + a, k, 1}) = phir[a].imag();
        out.Row({5 + a, k, 0}) = A[a].real();
        out.Row({5 + a, k, 1}) = A[a].imag();
        out.Row({9 + a, k, 0}) = V[a].real();
        out.Row({9 + a, k, 1}) = V[a].imag();
      }
    }
    out.fill();
  }

  out.close();
  in.close();
  
  H5Fclose(fileout_id);
  return 0;
}
