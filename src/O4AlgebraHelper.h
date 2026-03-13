#ifndef O4ALGEBRAHELPER
#define O4ALGEBRAHELPER

#include <petscdm.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <petscts.h>

#include <iostream>

#include <array>

// For the generator and other O4 conventions, see Derek's notes.

namespace O4AlgebraHelper {

PetscScalar norm(PetscScalar *v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void add(PetscScalar *a, PetscScalar *b, PetscScalar *c) {
  for (int i = 0; i < 3; ++i) {
    c[i] = a[i] + b[i];
  }
}

void scalmul(PetscScalar a, PetscScalar *b, PetscScalar *c) {
  for (int i = 0; i < 3; ++i) {
    c[i] = a * b[i];
  }
}

// c+=a*b
void muladd(PetscScalar a, PetscScalar *b, PetscScalar *c) {
  for (int i = 0; i < 3; ++i) {
    c[i] += a * b[i];
  }
}

void cross(PetscScalar *a, PetscScalar *b, PetscScalar *c) {
  c[0] = a[1] * b[2] - b[1] * a[2];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - b[0] * a[1];
}

PetscScalar dot(PetscScalar *a, PetscScalar *b) {
  PetscScalar res = 0;
  for (int i = 0; i < 3; ++i) {
    res += a[i] * b[i];
  }
  return res;
}

void quatmul(PetscScalar *v1, PetscScalar *v2, PetscScalar *w) {
  w[0] = v1[0] * v2[0] - dot(v1 + 1, v2 + 1);
  cross(v1 + 1, v2 + 1, w + 1);
  muladd(v2[0], v1 + 1, w + 1);
  muladd(v1[0], v2 + 1, w + 1);

  // res.getV() = v1[0] * v2.getV() + v2[0] * v1.getV() + cross(v1.getV(),
  // v2.getV());
}

void expmul(PetscScalar *v, PetscScalar *w) {
  PetscScalar abs = norm(v);
  if(abs > 1e-20){
    w[0] = cos(abs);
    scalmul(-(sin(abs) / abs), v, w + 1);
  }
  else{
    w[0] = 1;
    w[1] = 0.0;
    w[2] = 0.0;
    w[3] = 0.0;
  }
}

void O4Rotation(PetscScalar *V, PetscScalar *A,
                PetscScalar *phi) // In place O4 rotation specified by a vector
                                  // and axial charge
{
  PetscScalar tmp1[4], tmp2[4], tmp3[4];

  // Use the temporary to store L, R
  for (int i = 0; i < 3; ++i) {
    tmp1[i] = -0.5 * (V[i] - A[i]);
    tmp2[i] = 0.5 * (V[i] + A[i]);
  }

  // R rotation, in tmp3
  expmul(tmp1, tmp3);
  // L rotation in tmp1
  expmul(tmp2, tmp1);

  quatmul(phi, tmp3, tmp2);
  quatmul(tmp1, tmp2, phi);
}

}; // namespace O4AlgebraHelper

#endif
