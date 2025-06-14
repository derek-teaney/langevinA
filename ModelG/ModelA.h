#ifndef MODELASTRUCT
#define MODELASTRUCT

#include "NoiseGenerator.h"
#include "json/json.h"
#include <fstream>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <petscts.h>
#ifndef MODELA_NO_HDF5
#include <petscviewerhdf5.h>
#endif
#include "make_unique.h"

// This header file defines the structure and behavior of a simulation model
// called ModelA. It includes necessary libraries and defines several classes
// and structures to manage the simulation data, parameters, and operations.

// The ModelATime structure manages time-stepping information for the
// simulation. The ModelACoefficients structure holds thermodynamic and
// transport coefficients. The ModelAHandlerData structure contains options and
// settings for managing the simulation run. The ModelAData class aggregates
// all static data and configuration options for the simulation. The G_node and
// define data types for grid nodes in the simulation. The ModelA class
// encapsulates the entire simulation, including grid setup, initialization,
// and finalization.

// POD structure for recording the time stepping information. The finaltime is
// the final time of the simulation, the initialtime is the initial time, and
// deltat is the time step. The time variable is the current time.
struct ModelATime {
  PetscReal finaltime = 10.;
  PetscReal initialtime = 0.;
  PetscReal deltat = 0.24;
  PetscReal time = 0.;

  PetscReal t() const { return time; }
  PetscReal dt() const { return deltat; }
  PetscReal tinitial() const { return initialtime; }
  PetscReal tfinal() const { return finaltime; }
  void operator+=(const double &dtin) { time += dtin; }
  void reset() { time = initialtime; }

  // The finaltime is be adjusted so that the total number of time steps is an
  // integral number of saveFrequency*dt() units.  The new finaltime is
  // returned.
  PetscReal adjust_finaltime(const int &saveFrequency) {
    finaltime = initialtime + deltat * saveFrequency *
                                  static_cast<int>((finaltime - initialtime) /
                                                   (deltat * saveFrequency));
    return finaltime;
  }

  void read(const Json::Value &params) {
    finaltime = params.get("finaltime", finaltime).asDouble();
    initialtime = params.get("initialtime", initialtime).asDouble();
    deltat = params.get("deltat", deltat).asDouble();
    time = initialtime;
  }

  void print() {
    // Time Stepping
    PetscPrintf(PETSC_COMM_WORLD, "finaltime = %e\n", finaltime);
    PetscPrintf(PETSC_COMM_WORLD, "initialtime = %e\n", initialtime);
    PetscPrintf(PETSC_COMM_WORLD, "deltat  = %e\n", deltat);
  }
};

// POD structure for the thermodynamic and transport coefficients.
struct ModelACoefficients {

  PetscReal mass0 = -4.70052;
  PetscReal dmassdt = 0.;

  PetscReal lambda = 4.;
  PetscReal H = 0.003;
  PetscReal chi = 5;

  PetscReal gamma = 1.;
  PetscReal diffusion = 0.3333333;

  // Returns the value of the mass at a given time t.
  PetscReal mass(const double &t) const { return mass0 + dmassdt * t; }

  PetscReal sigma() const { return diffusion * chi; }
  PetscReal D() const { return diffusion; }

  void read(const Json::Value &params) {
    mass0 = params.get("mass0", mass0).asDouble();
    dmassdt = params.get("dmassdt", 0.).asDouble();

    lambda = params.get("lambda", lambda).asDouble();
    H = params.get("H", H).asDouble();
    chi = params.get("chi", chi).asDouble();

    gamma = params.get("gamma", gamma).asDouble();
    diffusion = params.get("diffusion", 1. / 3. * gamma).asDouble();
  }

  void print() {
    PetscPrintf(PETSC_COMM_WORLD, "mass0 = %e\n", mass0);
    PetscPrintf(PETSC_COMM_WORLD, "dmassdt = %e\n", dmassdt);

    PetscPrintf(PETSC_COMM_WORLD, "lambda = %e\n", lambda);
    PetscPrintf(PETSC_COMM_WORLD, "H = %e\n", H);
    PetscPrintf(PETSC_COMM_WORLD, "chi = %e\n", chi);
    PetscPrintf(PETSC_COMM_WORLD, "gamma = %e\n", gamma);
    PetscPrintf(PETSC_COMM_WORLD, "diffusion = %e\n", diffusion);
  }
};

// Lightweight data structure with access to options
struct ModelAHandlerData {

  std::string evolverType = "PV2HBSplit23";

  // random seed
  PetscInt seed = 10;

  // If we are to restore
  bool restart = false;

  // Options controlling the output. The outputfiletag labells the run. All
  // output files are tag_foo.txt, or tag_bar.h5
  std::string outputfiletag = "o4output";
  PetscInt saveFrequency = 3;
  PetscInt writeFrequency = -1;

  bool eventmode = false;
  int nevents = 1;
  int current_event = 0;
  double thermalization_time = 0.;

  // Quench Mode:
  //
  // In quench mode we start with some mass, quench_mode_mass0, and
  // thermalize the system with that mass. At time t=0 we start the simulation
  // with a different mass, as given by acoefficients.mass0
  bool quench_mode = false;
  double quench_mode_mass0 = -4.70052;

  // Initial amplitude, standing waves bool and dimension of init cond.
  PetscReal init_amp = 1.0;
  bool standing_waves = false;
  int init_dim = 1;

  void read(Json::Value &params) {
    evolverType = params.get("evolverType", evolverType).asString();
    seed = (PetscInt)params.get("seed", seed).asInt();
    restart = params.get("restart", false).asBool();
    outputfiletag = params.get("outputfiletag", "o4output").asString();
    saveFrequency = params.get("saveFrequency", saveFrequency).asInt();
    writeFrequency = params.get("writeFrequency", writeFrequency).asInt();
    thermalization_time =
        params.get("thermalization_time", thermalization_time).asDouble();

    quench_mode = params.get("quench_mode", quench_mode).asBool();
    quench_mode_mass0 =
        params.get("quench_mode_mass0", quench_mode_mass0).asDouble();

    eventmode = params.get("eventmode", eventmode).asBool();
    nevents = params.get("nevents", nevents).asInt();

    init_amp = params.get("init_amp", init_amp).asDouble();
    standing_waves = params.get("standing_waves", standing_waves).asBool();
    init_dim = params.get("init_dim", init_dim).asInt();
  }

  void print() {
    PetscPrintf(PETSC_COMM_WORLD, "evolverType = %s\n", evolverType.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "seed = %d\n", seed);
    PetscPrintf(PETSC_COMM_WORLD, "restart = %s\n",
                (restart ? "true" : "false"));
    PetscPrintf(PETSC_COMM_WORLD, "outputfiletag = %s\n",
                outputfiletag.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "saveFrequency = %d\n", saveFrequency);
    PetscPrintf(PETSC_COMM_WORLD, "thermalization_time = %e\n",
                thermalization_time);

    PetscPrintf(PETSC_COMM_WORLD, "quench_mode = %s\n",
                (quench_mode ? "true" : "false"));
    PetscPrintf(PETSC_COMM_WORLD, "quench_mode_mass0 = %e\n",
                quench_mode_mass0);

    PetscPrintf(PETSC_COMM_WORLD, "eventmode = %s\n",
                (eventmode ? "true" : "false"));
    PetscPrintf(PETSC_COMM_WORLD, "nevents = %d\n", nevents);
    PetscPrintf(PETSC_COMM_WORLD, "init_amp = %e\n", init_amp);
    PetscPrintf(PETSC_COMM_WORLD, "standing_waves = %s\n",
                (standing_waves ? "true" : "false"));
    PetscPrintf(PETSC_COMM_WORLD, "init_dim = %d\n", init_dim);
  }
};

// A lightweight data structure that contains all the information about the
// run. Every static bit of data should be accessible here
class ModelAData {

public:
  // Lattice dimension
  PetscInt NX = 16;
  PetscInt NY = 16;
  PetscInt NZ = 16;

  // Lattice size
  PetscReal LX = 16.;
  PetscReal LY = 16.;
  PetscReal LZ = 16.;

  // Lattice spacing in physical units
  PetscReal hX() const { return LX / NX; }
  PetscReal hY() const { return LY / NY; }
  PetscReal hZ() const { return LZ / NZ; }

  // Number of fields
  static const PetscInt Nphi = 4;
  static const PetscInt NA = 3;
  static const PetscInt NV = 3;
  static const PetscInt Ndof = Nphi + NA + NV;

  // Options controlling the clock
  ModelATime atime;

  // Thermodynamic and transport coefficients
  ModelACoefficients acoefficients;

  // Convenience function for returning the mass at the current time
  PetscReal mass() const { return acoefficients.mass(atime.t()); }

  // Options for management ;
  ModelAHandlerData ahandler;

public:
  ModelAData(Json::Value &params) {
    // Lattice. By default, NY=NX and NZ=NX.
    NX = params["NX"].asInt();
    NY = NX;
    NZ = NX;

    // By default, dx=dy=dz=1, namely LX=NX, LY=NY, LZ=NZ.
    LX = params.get("LX", NX).asInt();
    LY = LX;
    LZ = LX;

    // Time Stepping
    atime.read(params);
    acoefficients.read(params);
    ahandler.read(params);

    // Micro adjust the final time
    double finaltime = atime.adjust_finaltime(ahandler.saveFrequency);
    params["finaltime"] = finaltime;

    // Printout
    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) {
      print();
    }
  }

  //! Print out ModelAData for subsequent reading
  void print() {
    // Lattice
    PetscPrintf(PETSC_COMM_WORLD, "NX = %d\n", NX);
    PetscPrintf(PETSC_COMM_WORLD, "NY = %d\n", NY);
    PetscPrintf(PETSC_COMM_WORLD, "NZ = %d\n", NZ);
    PetscPrintf(PETSC_COMM_WORLD, "LX = %e\n", LX);
    PetscPrintf(PETSC_COMM_WORLD, "LY = %e\n", LY);
    PetscPrintf(PETSC_COMM_WORLD, "LZ = %e\n", LZ);

    // Time Stepping
    atime.print();

    // Transport and thermodynamic coefficients
    acoefficients.print();

    // Mangement parameters
    ahandler.print();
  }
};

/////////////////////////////////////////////////////////////////////////

// This describes the data at each lattice site of the grid
typedef struct {
  // this is the field phi for a = 1,2,3,4
  PetscScalar f[ModelAData::Nphi];
  // this is the the axial vector charge for s = 1,2,3
  PetscScalar A[ModelAData::NA];
  // this is the the vector charge for s = 1,2,3
  PetscScalar V[ModelAData::NV];
} G_node;

// This describes the data at each lattice site as a single vector x
typedef struct {
  PetscScalar x[ModelAData::Ndof];
} data_node;

// This is "the" class which contains access to all of the information
// about the run  as well as acess to the grid. It contains a copy of
// ModelAData which can be used to acess all configuration optons
class ModelA {

public:
  // Plain old data describing model A
  ModelAData data;

  // Domain descriptor
  DM domain;

  // Solution
  Vec solution;

  // Previous solution
  Vec previoussolution;

  // Rank of this processor
  int rank;

  // Construct the grid and initialize the fields according to  the
  // configuration parameters in the input ModelAData structure
  ModelA(const ModelAData &in) : data(in) {

    PetscInt stencil_width = 1;
    PetscInt Ndof = ModelAData::Nphi + ModelAData::NA + ModelAData::NV;
    DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                 DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR, data.NX, data.NY,
                 data.NZ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, Ndof,
                 stencil_width, NULL, NULL, NULL, &domain);

    DMSetFromOptions(domain);
    DMSetUp(domain);

    DMCreateGlobalVector(domain, &solution);
    VecDuplicate(solution, &previoussolution);

    // Setup the random number generation. If we are in in restart mode then we
    // try to read in the random number generator too.
    const auto &ahandler = data.ahandler;
    ModelARndm = make_unique<NoiseGenerator>(ahandler.seed);
    if (ahandler.restart and !ahandler.eventmode) {
      ModelARndm->read(ahandler.outputfiletag);
    }

    // Printout store the rank
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  }

  void finalize() {
    const auto &ahandler = data.ahandler;

    // When running in box mode we save the output and the
    // random number stream so we can restart the simulation
    // and build up statistics. In eventmode we are running
    // events, and there is no point in saving the state.
    if (not ahandler.eventmode) {
      ModelARndm->write(ahandler.outputfiletag);
      write(ahandler.outputfiletag);
    }
    VecDestroy(&previoussolution);
    VecDestroy(&solution);
    DMDestroy(&domain);
  }

  //! Reads in a stored initial condition from Derek's file.
  //! This is a helper function for initialize
  PetscErrorCode read(const std::string fnamein) {
    std::string fname = fnamein + "_save.h5";
#ifndef MODELA_NO_HDF5
    PetscViewer initViewer;
    PetscBool flg;
    PetscErrorCode ierr = PetscTestFile(fname.c_str(), '\0', &flg);
    if (!flg) {
      throw(
          std::string("Unable to open file in restart mode with filename = ") +
          fname);
    }
    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, fname.c_str(), FILE_MODE_READ,
                               &initViewer);
    CHKERRQ(ierr);
    PetscViewerSetFromOptions(initViewer);
    PetscObjectSetName((PetscObject)solution, "o4fields");
    ierr = VecLoad(solution, initViewer);
    CHKERRQ(ierr);
    PetscViewerDestroy(&initViewer);
    return ierr;
#else
    throw("ModelA::read Unable to load from file without HDF5 support. \n");
    return 0;
#endif
  }

  //! Reads in a stored initial condition from Derek's file.
  //! This is a helper function for initialize
  PetscErrorCode write(const std::string fnamein) {
    std::string fname = fnamein + "_save.h5";
#ifndef MODELA_NO_HDF5
    PetscViewer initViewer;
    PetscErrorCode ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, fname.c_str(),
                                              FILE_MODE_WRITE, &initViewer);
    CHKERRQ(ierr);
    PetscViewerSetFromOptions(initViewer);
    PetscObjectSetName((PetscObject)solution, "o4fields");
    ierr = VecView(solution, initViewer);
    CHKERRQ(ierr);
    PetscViewerDestroy(&initViewer);
    return 0;
#else
    PetscPrintf(
        "ModelA::write Unable to write to a file wihtout HDF5 support\n");
    return 0;
#endif
  }

  // Initialize the vector solution. If coldStart is false  then read in the
  // data from Derek's file. Otherwise fill with  random numbers, or if
  // zeroStart is true set to zero. Finally  if a function is provided,
  // f(x,y,z,L, params), this function will be used.
  PetscErrorCode initialize(double (*func)(const double &x, const double &y,
                                           const double &z, const int &L,
                                           void *params) = 0,
                            void *params = 0) {

    const auto &ahandler = data.ahandler;
    if (ahandler.restart) {
      try {
        read(ahandler.outputfiletag);
        return (0);
      } catch (const std::string &error) {
        std::cout << "Error in restart -- aborting!" << std::endl;
        std::cout << error << std::endl;
        std::abort();
      }
    }

    // Compute the lattice spacing
    PetscReal hx = data.hX();
    PetscReal hy = data.hY();
    PetscReal hz = data.hZ();

    // This Get a pointer to do the calculation
    PetscScalar ****u;
    PetscCall(DMDAVecGetArrayDOF(domain, solution, &u));

    // Get the Local Corner od the vector
    PetscInt i, j, k, L, xstart, ystart, zstart, xdimension, ydimension,
        zdimension;

    DMDAGetCorners(domain, &xstart, &ystart, &zstart, &xdimension, &ydimension,
                   &zdimension);

    // This is the actual computation of the thing
    for (k = zstart; k < zstart + zdimension; k++) {
      PetscReal z = k * hz;
      for (j = ystart; j < ystart + ydimension; j++) {
        PetscReal y = j * hy;
        for (i = xstart; i < xstart + xdimension; i++) {
          PetscReal x = i * hx;
          for (L = 0; L < ModelAData::Ndof; L++) {
            if (func) {
              u[k][j][i][L] = func(x, y, z, L, params);
            } else {
              u[k][j][i][L] = 0.;
            }
          }
        }
      }
    }
    PetscCall(DMDAVecRestoreArrayDOF(domain, solution, &u));
    return (0);
  }

  // Routine that initializes charge fields according to a gaussian distribution
  // and subtracts the total charge afterwards
  PetscErrorCode initialize_gaussian_charges() {
    // This Get a pointer to do the calculation
    PetscScalar ****u;
    PetscCall(DMDAVecGetArrayDOF(domain, solution, &u));

    // Get the Local Corner od the vector
    PetscInt i, j, k, L, xstart, ystart, zstart, xdimension, ydimension,
        zdimension;

    PetscCall(DMDAGetCorners(domain, &xstart, &ystart, &zstart, &xdimension,
                             &ydimension, &zdimension));

    // We are going initialize the grid with the charges being gaussian random
    // numbers. The charges are normalized so that the total charge is zero.
    std::vector<PetscScalar> charge_sum_local(ModelAData::Ndof, 0.);
    std::vector<PetscScalar> charge_sum(ModelAData::Ndof, 0.);

    PetscScalar chi = data.acoefficients.chi;
    for (k = zstart; k < zstart + zdimension; k++) {
      for (j = ystart; j < ystart + ydimension; j++) {
        for (i = xstart; i < xstart + xdimension; i++) {
          for (L = 0; L < ModelAData::Ndof; L++) {
            // Dont update the phi components
            if (L < ModelAData::Nphi) {
              continue;
            }

            // Generate gaussian random numbers for charges
            u[k][j][i][L] = sqrt(chi) * ModelARndm->normal();

            // Accumulate the total charge in a Buffer
            charge_sum_local[L] += u[k][j][i][L];
          }
        }
      }
    }

    // Find the total charge
    MPI_Allreduce(charge_sum_local.data(), charge_sum.data(), ModelAData::Ndof,
                  MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD);

    // Subtract the zero mode. Assumes lattice spacing is one
    PetscScalar V = data.NX * data.NY * data.NZ;
    for (k = zstart; k < zstart + zdimension; k++) {
      for (j = ystart; j < ystart + ydimension; j++) {
        for (i = xstart; i < xstart + xdimension; i++) {
          for (L = 0; L < ModelAData::Ndof; L++) {
            if (L < ModelAData::Nphi) {
              continue;
            }
            u[k][j][i][L] -= charge_sum[L] / V;
          }
        }
      }
    }

    PetscCall(DMDAVecRestoreArrayDOF(domain, solution, &u));

    return (0);
  }

  // Routine that initializes the fields randomly under the constraint phi^2 =
  // R, i.e. uniformly distributed spins on a 4d sphere
  //
  // Use hyperspherical coordinates:
  // s_0 = R * sin(phi) * sin(theta1) * sin(theta2)
  // s_1 = R * cos(phi) * sin(theta1) * sin(theta2)
  // s_2 = R * cos(theta1) * sin(theta2)
  // s_3 = R * cos(theta2)
  //
  // Measure: sin(theta1) sin^2(theta2) dphi dtheta1 dtheta2
  PetscErrorCode initialize_random_spins() {

    constexpr auto PI = 3.14159265358979323846;

    // This Get a pointer to do the calculation
    PetscScalar ****u;
    DMDAVecGetArrayDOF(domain, solution, &u);

    // Get the Local Corner od the vector
    PetscInt i, j, k, xstart, ystart, zstart, xdimension, ydimension,
        zdimension;

    DMDAGetCorners(domain, &xstart, &ystart, &zstart, &xdimension, &ydimension,
                   &zdimension);

    PetscReal R = data.ahandler.init_amp;
    PetscReal phi;
    PetscReal theta1;
    PetscReal theta2;

    // boolean and reals needed for rejection sampling
    bool accepted;
    PetscReal guess;
    PetscReal reference;

    // iterate over all lattice points
    for (k = zstart; k < zstart + zdimension; k++) {
      for (j = ystart; j < ystart + ydimension; j++) {
        for (i = xstart; i < xstart + xdimension; i++) {

          // get random uniform sample for phi
          phi = 2.0 * PI * ModelARndm->uniform();

          accepted = false;

          // loop for rejection sampling
          // (needed to sample thetas according to sin(theta1) * sin^2(theta2) )
          while (accepted == false) {
            // get random uniform sample for thetas and guess
            theta1 = PI * ModelARndm->uniform();
            theta2 = PI * ModelARndm->uniform();
            guess = ModelARndm->uniform();

            // calculate reference value
            reference = std::sin(theta1) * std::sin(theta2) * std::sin(theta2);

            // accept if guess <= reference, otherwise reject (loop again)
            if (guess <= reference)
              accepted = true;
          }

          // only initialize the field components
          u[k][j][i][0] =
              R * std::sin(phi) * std::sin(theta1) * std::sin(theta2);
          u[k][j][i][1] =
              R * std::cos(phi) * std::sin(theta1) * std::sin(theta2);
          u[k][j][i][2] = R * std::cos(theta1) * std::sin(theta2);
          u[k][j][i][3] = R * std::cos(theta2);
        }
      }
    }

    DMDAVecRestoreArrayDOF(domain, solution, &u);

    return (0);
  }

  // Routine that initializes the fields randomly according to a wave with wave
  // number k = 4pi/L and normalization phi^2 = init_amp = R
  //
  // if init_dim = 1 and standing_waves = false
  // s_0 = R * cos(kx)
  // s_1 = R * sin(kx)
  // n_01 = k * chi = omega * chi
  // rest 0
  //
  // if init_dim = 2
  // s_0 = R * cos(kx) cos(ky)
  // s_1 = R * sin(kx) cos(ky)
  // s_2 = R * sin(ky)
  // n_01 = k * chi = omega * chi
  // n_02 = k * chi = omega * chi
  // n_12 = k * chi = omega * chi
  // rest 0
  //
  // if init_dim = 1 and standing_waves = true
  // s_0 = R * cos(kx)
  // s_1 = 0
  // s_2 = R * sin(ky)
  // n_01 = k * chi = omega * chi
  // rest 0
  PetscErrorCode initialize_wave_spins() {

    constexpr auto PI = 3.14159265358979323846;

    PetscScalar chi = data.acoefficients.chi;
    PetscReal wave_k = 4 * PI / data.LX;
    PetscReal argument;

    // This Get a pointer to do the calculation
    PetscScalar ****u;
    DMDAVecGetArrayDOF(domain, solution, &u);

    // We are going initialize the grid with the charges being gaussian random
    // numbers. The charges are normalized so that the total charge is zero.
    std::vector<PetscScalar> charge_sum_local(ModelAData::Ndof, 0.);
    std::vector<PetscScalar> charge_sum(ModelAData::Ndof, 0.);

    // Get the Local Corner od the vector
    PetscInt i, j, k, L, xstart, ystart, zstart, xdimension, ydimension,
        zdimension;

    DMDAGetCorners(domain, &xstart, &ystart, &zstart, &xdimension, &ydimension,
                   &zdimension);
    for (k = zstart; k < zstart + zdimension; k++) {
      for (j = ystart; j < ystart + ydimension; j++) {
        for (i = xstart; i < xstart + xdimension; i++) {

          argument = wave_k;
          for (L = 0; L < ModelAData::Ndof; L++) {

            // field components s_a
            if (L < ModelAData::Nphi) {

              // s_0 component
              if (L == 0) {
                u[k][j][i][L] = data.ahandler.init_amp;
                // 1d init cond.
                if (data.ahandler.init_dim == 1) {
                  u[k][j][i][L] *= cos(argument * i);
                }
                // 2d init cond.
                else if (data.ahandler.init_dim == 2) {
                  u[k][j][i][L] *= cos(argument * i) * cos(argument * j);
                }
              }
              // s_1 component
              else if (L == 1) {
                u[k][j][i][L] = data.ahandler.init_amp;
                // 1d init cond.
                if (data.ahandler.init_dim == 1) {
                  // if standing wave solution
                  if (data.ahandler.standing_waves) {
                    u[k][j][i][L] *= 0.0;
                  } else {
                    u[k][j][i][L] *= sin(argument * i);
                  }
                }
                // 2d init cond.
                else if (data.ahandler.init_dim == 2) {
                  u[k][j][i][L] *= sin(argument * i) * cos(argument * j);
                }
              }
              // s_2 component
              else if (L == 2) {
                u[k][j][i][L] = data.ahandler.init_amp;
                // 1d init cond.
                if (data.ahandler.init_dim == 1) {
                  // if standing wave solution
                  if (data.ahandler.standing_waves) {
                    u[k][j][i][L] *= sin(argument * i);
                  } else {
                    u[k][j][i][L] *= 0.0;
                  }
                }
                // 2d init cond.
                else if (data.ahandler.init_dim == 2) {
                  u[k][j][i][L] *= sin(argument * j);
                }
              }
            }

            // charges n_A (4,5,6) and n_V (7,8,9)
            else {
              // (n_A)_0 (or n_01)
              if (L == 4) {
                u[k][j][i][L] = wave_k * chi;
              }
              // (n_A)_0 (or n_01) OR (n_V)_2 (or n_12)
              else if (L == 5 || L == 9) {
                // 2d init cond.
                if (data.ahandler.init_dim == 2) {
                  u[k][j][i][L] = wave_k * chi;
                }
              }
            }
          } // end of assignment for single point
        }
      }
    } // end of loop over all points

    DMDAVecRestoreArrayDOF(domain, solution, &u);

    return (0);
  }
};
#endif
