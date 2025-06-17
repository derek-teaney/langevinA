#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "ModelA.h"
#include "NoiseGenerator.h"
#include "Stepper.h"
#include "gitversion.h"
#include "make_unique.h"

// Measurer, where the Petsc are included
#include "measurer.h"
#include "measurer_output.h"

void thermalize_event(ModelA *const model) {
  const auto &ahandler = model->data.ahandler;
  auto &atime = model->data.atime;
  auto &acoefficients = model->data.acoefficients;

  // Initialize a quench.  Set the initial temperature (mass parameter)
  // according a given value and thermalize this initial condition.  Then,
  // after the thremalization process reset the mass to the one used for the
  // actual running (as opposed to initializing) the code. The reset process
  // is handled below
  const double mass0 = acoefficients.mass0; // Store the mass for reset process
  const double dmassdt =
      acoefficients.dmassdt; // Store the slope for reset process
  if (ahandler.quench_mode) {
    acoefficients.mass0 = ahandler.quench_mode_mass0;
    acoefficients.dmassdt = 0.;
    PetscPrintf(PETSC_COMM_WORLD,
                "Settinng up a quench initial condition with initial mass %e\n",
                acoefficients.mass0);
  }

  // Thermalize the state in memory at the initial time ;
  int nsteps = static_cast<int>(ahandler.thermalization_time / atime.dt());
  PetscPrintf(PETSC_COMM_WORLD, "Thermalizing event %d\n",
              ahandler.current_event);

  // Thermalize the initialconditions
  std::unique_ptr<EulerLangevinHB> thermalizer =
      std::make_unique<EulerLangevinHB>(*model);
  model->initialize_gaussian_charges();
  for (int i = 0; i < nsteps; i++) {
    const int substeps = 6;
    for (int j = 0; j < substeps; j++) {
      thermalizer->step(atime.dt() / substeps);
    }
    PetscPrintf(PETSC_COMM_WORLD,
                "Thermalizing Event/Timestep %d/%d: step size = %g, time = %g, "
                "nsteps to thermalize = %d, mass %e \n",
                ahandler.current_event, i, (double)atime.dt(),
                (double)atime.t(), nsteps, model->data.mass());
  }
  thermalizer->finalize();

  // If we are performing a quench, set the mass back to its nominal value.
  if (ahandler.quench_mode) {
    PetscPrintf(PETSC_COMM_WORLD,
                "Finalizing  quench initial condition with initial mass %e\n",
                acoefficients.mass0);

    acoefficients.mass0 = mass0;
    acoefficients.dmassdt = dmassdt;

    PetscPrintf(PETSC_COMM_WORLD, "and final initial mass %e\n",
                acoefficients.mass0);
  }
}

void run_event(ModelA *const model, Stepper *const step) {

  const auto &ahandler = model->data.ahandler;
  auto &atime = model->data.atime;
  atime.reset();

  if (not ahandler.restart) {
    thermalize_event(model);
  }

  // Set up logging for PETSc so we can find out how much time
  // each part takes
  PetscInt steps = 0;
  PetscLogEvent measurements, stepmonitor, saving;
  PetscLogEventRegister("Measurements", 0, &measurements);
  PetscLogEventRegister("Saving the fields", 0, &saving);
  PetscLogEventRegister("Steps", 0, &stepmonitor);

  // Set filename
  std::string filename;
  if (ahandler.eventmode) {
    std::stringstream namestream;
    namestream << ahandler.outputfiletag << "_" << std::setw(4)
               << std::setfill('0') << ahandler.current_event << ".h5";
    filename = namestream.str();
  } else {
    filename = ahandler.outputfiletag + ".h5";
  }
  // Set file access
  PetscFileMode file_access = FILE_MODE_WRITE;
  if (ahandler.restart) {
    file_access = FILE_MODE_APPEND;
  }
  // Open the file and create the measurement object
  Measurer measurer(model);

  int rank = -1;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  std::unique_ptr<measurer_output_fasthdf5> measurer_output;
  if (rank == 0) {
    measurer_output = std::make_unique<measurer_output_fasthdf5>(
        &measurer, filename, file_access);
  }

  // Start the loop
  const double tiny = 1.e-10;
  while (atime.t() < atime.tfinal() - tiny) {

    // measure the solution every saveFrequency
    if (steps % ahandler.saveFrequency == 0) {
      PetscLogEventBegin(measurements, 0, 0, 0, 0);
      measurer.measure(&model->solution);
      if (rank == 0) {
        measurer_output->save();
      }
      PetscPrintf(PETSC_COMM_WORLD,
                  "Event/Timestep %d/%d: step size = %g, time = %g, final = "
                  "%g, mass = %e\n",
                  ahandler.current_event, steps, (double)atime.dt(),
                  (double)atime.t(), (double)atime.tfinal(),
                  model->data.mass());
      PetscLogEventEnd(measurements, 0, 0, 0, 0);
    }

    // Write the solution to tape if writeFrequency > 0. This is used for
    // plotting of the solution. It is normally not analyzed, or written.
    if (ahandler.writeFrequency > 0 and steps % ahandler.writeFrequency == 0) {
      PetscLogEventBegin(saving, 0, 0, 0, 0);
      std::ostringstream tString;
      tString << std::setprecision(4) << "_t_" << atime.t();
      model->write(ahandler.outputfiletag + tString.str());
      PetscLogEventEnd(saving, 0, 0, 0, 0);
    }

    // Do the actual steps
    PetscLogEventBegin(stepmonitor, 0, 0, 0, 0);
    step->step(atime.dt());
    PetscLogEventEnd(stepmonitor, 0, 0, 0, 0);

    // Increment the clock
    steps++;
    atime += atime.dt();
  }
}

void Run(nlohmann::json &inputs) {
  // Digest the inputs, some fields may be modified on ouptut
  ModelAData inputdata(inputs);

  // allocate the grid and initialize
  ModelA model(inputdata);

  // Read in the initial conditions or initialize to zero
  model.initialize();

  // Construct the stepper
  std::unique_ptr<Stepper> step;
  auto &etype = inputdata.ahandler.evolverType;
  if (etype == "PV2HBSplit23") {
    std::array<unsigned int, 2> s = {2, 3};
    // Default is to include all steps
    step = std::make_unique<PV2HBSplit>(model, s);
  } else if (etype == "PV2HBSplit23NoDiffuse") {
    std::array<unsigned int, 2> s = {2, 3};
    const bool ideal = true;
    const bool heatbath = true;
    const bool diffuse = false;
    step = std::make_unique<PV2HBSplit>(model, s, ideal, heatbath, diffuse);
  } else if (etype == "PV2HBSplit23OnlyDiffuse") {
    std::array<unsigned int, 2> s = {2, 3};
    const bool ideal = false;
    const bool heatbath = false;
    const bool diffuse = true;
    step = std::make_unique<PV2HBSplit>(model, s, ideal, heatbath, diffuse);
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "Unrecognized stepper type %s. Aborting...\n",
                etype.c_str());
    return ;
  }

  auto &ahandler = model.data.ahandler;
  if (ahandler.eventmode) {
    for (int i = 0; i < ahandler.nevents; i++) {
      run_event(&model, step.get());
      ahandler.current_event++;
    }
  } else {
    run_event(&model, step.get());
  }
  // Destroy everything
  step->finalize();
  model.finalize();
}


int main(int argc, char **argv) {

  // Initialization of PETSc universe
  PetscErrorCode ierr;
  std::string help =
      "Usage: \n\n\t ./SuperPions.exe -input input.json [options]\n\n";
  ierr = PetscInitialize(&argc, &argv, (char *)0, help.c_str());
  if (ierr) {
    return ierr;
  }

  // Open the input file and parse the inputs into the ModelAData
  char filename[PETSC_MAX_PATH_LEN] = "";
  ierr = PetscOptionsGetString(NULL, NULL, "-input", filename, sizeof(filename),
                               NULL);
  nlohmann::json  inputs;
  std::ifstream ifs(filename);
  if (ifs) {
    ifs >> inputs;
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "Unable to open input file %s. Aborting...\n",
                filename);
    return PetscFinalize();
  }
  PetscPrintf(PETSC_COMM_WORLD, "Current version: %s\n", gitversion);
  
  Run(inputs);

  return PetscFinalize();
}
