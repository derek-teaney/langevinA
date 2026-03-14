#include <iostream>
#include "O4AlgebraHelper.h"

void print_rotation(
      std::array<PetscScalar,3> &V,
      std::array<PetscScalar,3> &A,
      std::array<PetscScalar,4> &phi) 
{
  O4AlgebraHelper::O4Rotation(V.data(), A.data(), phi.data())  ;
  for (int i = 0; i < 4; ++i) {
    std::cout << phi[i] << std::endl;
  }
}


// This is a test of the sign conventions of the O4AlgebraHelper.
// All of these steps should be consistent with the SO(4) note.
int main(int argc, char **argv) 
{
  std::array<PetscScalar,3>A{M_PI/4., 0, 0};
  std::array<PetscScalar,3>V{0, 0, 0};
  std::array<PetscScalar,4>phi{1, 0, 0, 0};
  std::cout << "Axial of (1, 0, 0, 0) rotation by  Pi/4 theta_{Ax}" << std::endl ;
  std::cout << "Our conventions are this should be  (1/sqrt{2}, -1/sqrt{2}, 0, 0)" << std::endl ;
  print_rotation(V, A, phi) ;

  A = {0, 0, 0} ;
  V = {0, 0, M_PI/4} ;
  phi = {0, 1, 0, 0} ;
  std::cout << "Vector of (0, 1, 0, 0) rotation by  Pi/4 theta_{Vz}" << std::endl ;
  std::cout << "Our conventions are this should be  (0, 1/sqrt{2}, -1/sqrt{2}, 0)" << std::endl ;
  print_rotation(V, A, phi) ;

  A = {0, 0, 0} ;
  V = {0, 0, M_PI/4} ;
  phi = {0, 0, 1, 0} ;
  std::cout << "Vector of (0, 0, 1, 0) rotation by  Pi/4 theta_{Vz}" << std::endl ;
  std::cout << "Our conventions are this should be  (0, 1/sqrt{2}, 1/sqrt{2}, 0)" << std::endl ;
  print_rotation(V, A, phi) ;

}
