#include "NoiseGenerator.h"
#include <petscsys.h>

int main(int argc, char **argv) {

  PetscErrorCode ierr;
  std::string help = "Random number tests";
  ierr = PetscInitialize(&argc, &argv, (char *)0, help.c_str());
  if (ierr) {
    return ierr;
  }
  int rank, size;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  // Generate and  Write
  int baseSeed = 2021;
  NoiseGenerator engine_w(baseSeed);
  PetscPrintf(PETSC_COMM_WORLD, "Initializing the engine with a baseseed and "
                                "test the generation on multiprocessors\n");
  PetscPrintf(PETSC_COMM_WORLD, "[counter][rank] number\n");
  for (size_t i = 0; i < 4; i++) {
    PetscScalar x = engine_w.uniform();
    if (i == 1) {
      engine_w.write("test_random");
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d][%d] %f\n", i, rank, x);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }

  // Read and generate
  PetscPrintf(PETSC_COMM_WORLD,
              "Read in the random number state and generate\n");
  PetscPrintf(PETSC_COMM_WORLD, "[counter][rank] number\n");
  NoiseGenerator engine_r;
  engine_r.read("test_random");
  for (size_t i = 2; i < 4; i++) {
    PetscScalar x = engine_r.uniform();
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "x[%d][%d] %f\n", i, rank, x);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }
  // So the printing is not screwed up
  MPI_Barrier(PETSC_COMM_WORLD);

  NoiseGenerator engine(baseSeed + 1);

  PetscPrintf(
      PETSC_COMM_WORLD,
      "Clocking the uniform random number generator... use -log_view\n");
  PetscLogEvent uniformGeneration;
  const size_t N = 1000000000;
  PetscLogEventRegister("Uniform", 0, &uniformGeneration);
  // Tests speed uniform
  PetscLogEventBegin(uniformGeneration, 0, 0, 0, 0);
  double xx = 0.;
  for (size_t i = 0; i < N; i++) {
    PetscScalar x = engine.variance1();
    xx += x * x;
    if (i < 2) {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d][%d] %f\n", i, rank, x);
      PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    }
  }
  xx /= static_cast<double>(N);
  PetscLogEventEnd(uniformGeneration, 0, 0, 0, 0);
  PetscPrintf(PETSC_COMM_WORLD, "The variance is %g\n", xx);

  // So the printing is not screwed up
  MPI_Barrier(PETSC_COMM_WORLD);

  // Test speed normal
  PetscPrintf(
      PETSC_COMM_WORLD,
      "Clocking the gaussian random number generator... use -log_view\n");
  PetscLogEvent normalGeneration;
  PetscLogEventRegister("Normal", 0, &normalGeneration);
  PetscLogEventBegin(normalGeneration, 0, 0, 0, 0);
  xx = 0.;
  for (size_t i = 0; i < N; i++) {
    PetscScalar x = engine.normal();
    xx += x * x;
    if (i < 2) {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d][%d] %f\n", i, rank, x);
      PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    }
  }
  xx /= static_cast<double>(N);
  PetscLogEventEnd(normalGeneration, 0, 0, 0, 0);
  PetscPrintf(PETSC_COMM_WORLD, "The variance is %g\n", xx);

  PetscFinalize();

  return 0;
}
