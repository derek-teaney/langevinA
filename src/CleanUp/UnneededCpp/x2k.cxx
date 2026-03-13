
#include "hdf5.h"
#include "ntuple.h"
#include <array>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <vector>

inline bool exists_test(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

// The purpose of the program is to take the fourier transform of a "wall"
// like dataset.  Example:
//
// ./x2k.exe foo.h5 wally
//
// This will open foo.h5, and create foo_out.h5 if does not exist.  And
// then in foo_out.h5 will create a dataset,  wally_k which is a fourier
// transform of the original data set.
//
// wally is assumed to be of the form:
//
// wx = wally[ievent, jfield, x]
//
// Where the first index is the event number, the second index labels the
// field, and x labels the lattice site, with the total number of sites an
// even number. The output takes the form
//
// wk = wally_k[ievent, jfield, k, 0 and 1] = 1/N \sum_{k} e^{-i 2pi k x/N} wx
//
// The complex numbers are stored sequetially  with
// 0/1 being the the real/imag parts
//
// This is accessed from numpy and hdf5. For example to
// retrieve a column of complex numbers do the following
//
// data = filehhandle["wally_k""][:,0,1,:]
//
// cdata = data.view(dtype=np.complex128)
int main(int argc, char **argv) {

  if (argc < 3) {
    std::cout << "Usage x2k filename.h5 datasetname" << std::endl;
    exit(0);
  }

  std::string fin(argv[1]);
  std::string dsetin(argv[2]);

  std::string fout = fin;
  fout = fout.erase(fout.find_last_of("."), std::string::npos) + "_out.h5";

  std::string dsetout = dsetin + "_k";
  std::cout << "Creating dataset " << dsetout << " in " << fout << std::endl;

  // Open the filespace and grab the ntuples
  hid_t filein_id = H5Fopen(fin.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  ntuple<2> in(dsetin, filein_id);

  // get the dimensions of the input tuple
  auto N = in.getN();

  // open the output file
  hid_t fileout_id;
  if (exists_test(fout)) {
    fileout_id = H5Fopen(fout.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  } else {
    fileout_id =
        H5Fcreate(fout.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }

  // Remove the old output if it exists
  if (H5Lexists(fileout_id, dsetout.c_str(), H5P_DEFAULT)) {
    H5Ldelete(fileout_id, dsetout.c_str(), H5P_DEFAULT);
  }

  // Construct the dimensions of the output tuple
  std::array<size_t, 3> NK{N[0], N[1] / 2 + 1, 2};

  ntuple<3> out(NK, dsetout, fileout_id);

  // Create the storage space for the plan
  std::vector<double> in_ptr(N[1], 0.);
  std::vector<std::complex<double>> out_ptr(NK[1], 0.);

  // Create the plane for the fourier transform
  fftw_plan p = fftw_plan_dft_r2c_1d(
      N[1], in_ptr.data(), reinterpret_cast<fftw_complex *>(out_ptr.data()),
      FFTW_MEASURE);

  for (size_t i = 0; i < in.nrows(); i++) {
    in.readrow(i);
    for (size_t i0 = 0; i0 < N[0]; i0++) {
      // Copy the data to the in_ptr of the fouier transform
      for (size_t i1 = 0; i1 < N[1]; i1++) {
        size_t k = in.at({i0, i1});
        in_ptr[i1] = in.row[k];
      }

      fftw_execute(p);

      // Copy the data to the out_ptr of the fourier transform
      for (size_t i1 = 0; i1 < NK[1]; i1++) {
        size_t k = out.at({i0, i1, 0});
        out.row[k] = out_ptr[i1].real()/N[1];
        out.row[k + 1] = out_ptr[i1].imag()/N[1];
      }
    }
    out.fill();
  }
  fftw_destroy_plan(p);

  // clean up
  out.close();
  in.close();
  H5Fclose(filein_id);
  H5Fclose(fileout_id);
  return 0;
}
