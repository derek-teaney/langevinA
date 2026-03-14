#include "hdf5.h"
#include "ntuple.h"
#include <array>

int main(int argc, char **argv) {

  hid_t file_id;
  file_id = H5Fcreate("work.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  static const hsize_t N = 8;
  std::array<double, N> datain = {1, 2, 3, 4, 5, 6, 7, 8};

  std::array<size_t, 1> NN1 = {8};
  ntuple<1> nt1(NN1, "phi1", file_id);

  std::array<size_t, 2> NN2 = {2, 4};
  ntuple<2> nt2(NN2, "phi2", file_id);

  std::array<size_t, 3> NN3 = {2, 2, 2};
  ntuple<3> nt3(NN3, "phi3", file_id);

  for (size_t j = 0; j < 3; j++)  {

    for (size_t i0 = 0; i0 < NN1[0]; i0++) {
      // Seems a bit excessive in 1d
      // size_t k = nt.at({i0, i1}) ;
      nt1.row[i0] = datain[i0];
    }

    for (size_t i0 = 0; i0 < NN2[0]; i0++) {
      for (size_t i1 = 0; i1 < NN2[1]; i1++) {
        size_t k = nt2.at({i0, i1});
        nt2.row[k] = datain[k];
      }
    }

    for (size_t i0 = 0; i0 < NN3[0]; i0++) {
      for (size_t i1 = 0; i1 < NN3[1]; i1++) {
        for (size_t i2 = 0; i2 < NN3[2]; i2++) {
          size_t k = nt3.at({i0, i1, i2});
          nt3.row[k] = datain[k];
        }
      }
    }
    nt1.fill();
    nt2.fill();
    nt3.fill();
  }

  // Close the ntuple
  nt1.close();
  nt2.close();
  nt3.close();

  // close the filespace
  H5Fclose(file_id);

  // Open the filespace and grab the ntuples
  file_id = H5Fopen("work.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  ntuple<1> in1("phi1", file_id);
  ntuple<2> in2("phi2", file_id);
  ntuple<3> in3("phi3", file_id);

  // Get the dimensioning
  auto IN1 = in1.getN();
  auto IN2 = in2.getN();
  auto IN3 = in3.getN();

  auto D2 = in2.get_dims();
  printf("The dimensions are: ");
  for (size_t i = 0; i < 3; i++) {
    printf("%lu ", (D2[i]));
  }
  printf("\n");

  printf("Ntumple dimensions are : ");
  for (size_t i = 0; i < 2; i++) {
    printf("%lu ", IN2[i]);
  }
  printf("\n");

  for (size_t i = 0; i < in1.nrows(); i++) {
    in1.readrow(i);
    for (size_t j = 0; j < IN1[0]; j++) {
      printf("%15.5e ", in1.row[j]);
    }
    printf("\n");

    printf("\n");
    in2.readrow(i);
    for (size_t i0 = 0; i0 < IN2[0]; i0++) {
      for (size_t i1 = 0; i1 < IN2[1]; i1++) {
        size_t k = nt2.at({i0, i1});
        printf("%15.5e ", in2.row[k]);
      }
      printf("\n");
    }
    printf("\n");

    printf("\n");
    in3.readrow(i);
    for (size_t i0 = 0; i0 < IN3[0]; i0++) {
      for (size_t i1 = 0; i1 < IN3[1]; i1++) {
        for (size_t i2 = 0; i2 < IN3[2]; i2++) {
          size_t k = nt3.at({i0, i1, i2});
          printf("%15.5e ", in3.row[k]);
        }
        printf("\n");
      }
      printf("\n");
    }
    printf("\n");
  }

  in1.close();
  in2.close();
  in3.close();

  // close the filespace
  H5Fclose(file_id);
  return 0;
}
