#ifndef NOISEGENERATOR
#define NOISEGENERATOR

#include <fstream>
#include <iostream>
#include <memory>
#include <petscsys.h>
#include <petscvec.h>
#include <random>
#include <string>
#include <vector>

#define NOISEGENERATOR_STDCPP

#ifdef NOISEGENERATOR_STDCPP
#include "xoroshiro128plus.h"
#endif

#ifdef NOISEGENERATOR_GSL
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

class gsl_rng_urng {
public:
  typedef unsigned long int result_type;
  gsl_rng *grng;

  static constexpr unsigned long int min() { return 0; }
  static constexpr unsigned long int max() { return 0xFFFFFFFF; }
  gsl_rng_urng() { grng = gsl_rng_alloc(gsl_rng_mt19937); }
  ~gsl_rng_urng() { gsl_rng_free(grng); }
  void seed(unsigned long int sd) { gsl_rng_set(grng, sd); }

  unsigned long int operator()() { return gsl_rng_get(grng); }
  void write(const std::string &fname) {
    FILE *fp = fopen(fname.c_str(), "wb");
    if (!fp) {
      std::cout << "gsl_rng_urng::write Unable to open file for writing in gsl URNG." << std::endl;
    } else {
      gsl_rng_fwrite(fp, grng);
    }
    fclose(fp);
  }
  void read(const std::string &fname) {
    FILE *fp = fopen(fname.c_str(), "rb");
    if (!fp) {
      std::cout << "gsl_rng_urng::read Unable to open file for reading in gsl URNG. Defaulting to seed." << std::endl;
      return ;
    } else {
      gsl_rng_fread(fp, grng);
    }
    fclose(fp);
  }
};
#endif

// typedef std::ranlux48 RNGType;
// typedef std::mt19937_64 RNGType;
typedef xoroshiro128plus RNGType;
// typedef gsl_rng_urng RNGType;

class NoiseGenerator {
public:
  // Initialize the random number generator according to a base seed.
  //
  // Each processor receives its own seed (based on the base seed), uisng
  // the c++ seed_seq function
  NoiseGenerator(const int &baseSeed = 0)
      : uniformDistributionRt3(-sqrt(3.), sqrt(3.)) {

    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int size = 0;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    std::seed_seq sequence = {baseSeed};
    std::vector<unsigned int> seeds(size, 0);
    sequence.generate(seeds.begin(), seeds.end());
    // Take the seed according to our rank
    rng.seed(seeds[rank]);
  }

  ~NoiseGenerator() {}

  PetscErrorCode fillVec(Vec U) {
    PetscScalar *array;
    PetscInt nloc = 0;
    VecGetLocalSize(U, &nloc);
    VecGetArray(U, &array);

    for (int i = 0; i < nloc; i++) {
      array[i] = normal();
    }

    VecRestoreArray(U, &array);
    return 0;
  }

  PetscReal normal() { return normalDistribution(rng); }

  PetscReal uniform() { return uniformDistribution(rng); }

  PetscReal uniformRt3() { return uniformDistributionRt3(rng); }

  PetscReal variance1() { return uniformDistributionRt3(rng); }

  RNGType &generator() { return rng; }

  //! Write out the state of the random number generator
  //! for each processor, the output files take the form
  //! name_0.rng, name_1.rng, ...
  void write(const std::string &filename_stub) {
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    std::string srank = std::to_string(rank);
    std::string fname = filename_stub + "_" + srank + ".rng";
#ifdef NOISEGENERATOR_GSL
    rng.write(fname);
#else
    std::ofstream out(fname);
    if (out) {
      out << rng;
    } else {
      std::cout << "urng::write Unable to open file for writing in URNG. Defaulting to seed." << std::endl;
    }
#endif
  }

  // Read back in the random number generators written by  write. If the file
  //does not exist the default seed  is used. See the constructor
  void read(const std::string &filename_stub) {
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    std::string srank = std::to_string(rank);
    std::string fname = filename_stub + "_" + srank + ".rng";
#ifdef NOISEGENERATOR_GSL
    rng.read(fname);
#else
    std::ifstream in(fname);
    if (in) {
      in >> rng;
    } else {
      std::cout << "urng::read Unable to open file for reading in URNG. Defaulting to seed." << std::endl;
    }
#endif
  }

private:
  RNGType rng;
  // Normal distribution with variance 1
  std::normal_distribution<PetscReal> normalDistribution;
  // Flat distribution between -sqrt(3) .. sqrt(3) with variance 1
  std::uniform_real_distribution<PetscReal> uniformDistributionRt3;
  // Flat distribution between 0 and 1
  std::uniform_real_distribution<PetscReal> uniformDistribution;
};

//! Global random number generator which is initialized when
//! ModelA is created.
extern std::unique_ptr<NoiseGenerator> ModelARndm;
#endif
