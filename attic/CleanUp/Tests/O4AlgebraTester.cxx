#include "O4AlgebraHelper.h"

#include <iostream>

int main() {

  PetscScalar a[3]{1, 2, 3};
  PetscScalar b[3]{4, 6, 7};
  PetscScalar c[3], d[3], q[4];

  O4AlgebraHelper::add(a, b, c);
  O4AlgebraHelper::scalmul(2, b, d);

  for (int i = 0; i < 3; ++i) {
    std::cout << c[i] << std::endl;
    std::cout << d[i] << std::endl;
  }

  O4AlgebraHelper::expmul(a, q);

  for (int i = 0; i < 4; ++i) {
    std::cout << q[i] << std::endl;
  }

  PetscScalar v1[4]{0.9, 0.6, 0.3, 0.5}, v2[4]{0.1, 0.2, 0.6, 0.9}, v3[4];

  O4AlgebraHelper::quatmul(v1, v2, v3);

  for (int i = 0; i < 4; ++i) {
    std::cout << v3[i] << std::endl;
  }

  PetscScalar A[3]{1, 2, 3};
  PetscScalar V[3]{4, 5, 6};
  PetscScalar phi[4]{0.5, 0.5, 0.5, 0.5};

  O4AlgebraHelper::O4Rotation(V, A, phi);

  for (int i = 0; i < 4; ++i) {
    std::cout << phi[i] << std::endl;
  }

  return 0;
}
