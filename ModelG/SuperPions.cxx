#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

#include "ModelA.h"
#include "Stepper.h"
#include "gitversion.h"
#include "initialize.h"

// Measurer, where the Petsc are included
#include "measurer.h"
#include "measurer_output.h"

void thermalize_event(ModelA *const model) {
  const auto &ahandler = model->data.ahandler;
  auto &atime = model->data.atime;

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
}

void initialize_event(const int &ievent, ModelA *const model,
                      nlohmann::json &inputs) {
  const auto &ahandler = model->data.ahandler;
  std::string initialization = inputs["initialization"];
  bool zero_charges = inputs.value("zero_charges", false);
  
  if (initialization == "default") {
    // Do a cold start and thermalize the event
    if (ievent == 0) {
      model->initialize();
    }
    thermalize_event(model);

  } else if (initialization == "restart") {
    // Look for a previously saved initial condtitions
    model->read(ahandler.outputfiletag);

  } else if (initialization == "quench_mode") {
    // Initialize a quench.  Set the initial temperature (mass parameter)
    // according a given value and thermalize this initial condition.  Then,
    // after the thremalization process reset the mass to the one used for the
    // actual running (as opposed to initializing) the code. The reset process
    // is handled below

    auto &acoefficients = model->data.acoefficients;
    const double mass0 =
        acoefficients.mass0; // Store the mass for reset process
    const double dmassdt =
        acoefficients.dmassdt; // Store the slope for reset process

    // Set the quench mass
    acoefficients.mass0 = ahandler.quench_mode_mass0;
    acoefficients.dmassdt = 0.;
    PetscPrintf(PETSC_COMM_WORLD,
                "Setting up a quench initial condition with initial mass %e\n",
                acoefficients.mass0);

    // Thermalize at the quench mass
    if (ievent == 0) {
      model->initialize();
    }
    thermalize_event(model);

    // Reset the mass and teh slope
    acoefficients.mass0 = mass0;
    acoefficients.dmassdt = dmassdt;

    PetscPrintf(PETSC_COMM_WORLD, "and final initial mass %e\n",
                acoefficients.mass0);

  } else if (initialization == "randomspins") {
    model->initialize_random_spins();
    if (zero_charges) {
        PetscScalar zero = 0.0;
        model->initialize_gaussian_const(&zero);
    }
    else {
        model->initialize_gaussian_charges();
    }

  } else if (initialization == "gaussians") {
    model->initialize(initialize_gaussians, &inputs["gaussians"]);
    model->write(inputs["outputfiletag"].get<std::string>() + "_initial");

  } else if (initialization == "spinwaves") {
    model->initialize(initialize_wave_spins, &inputs["spinwaves"]);
    model->write(inputs["outputfiletag"].get<std::string>() + "_initial");

  } else {
    throw std::runtime_error(
        "Unknown initialization type: " + initialization +
        ". Please use one of the following: default, restart, quench_mode, "
        "randomspins, gaussians, spinwaves.");
  }
}

void run_event(const int &ievent, ModelA *const model, Stepper *const step,
               nlohmann::json &inputs) {

  const auto &ahandler = model->data.ahandler;
  auto &atime = model->data.atime;

  // Set up logging for PETSc so we can find out how much time
  // each part takes
  int steps = 0;
  PetscLogEvent measurements, stepmonitor, saving;
  PetscLogEventRegister("Measurements", 0, &measurements);
  PetscLogEventRegister("Saving the fields", 0, &saving);
  PetscLogEventRegister("Steps", 0, &stepmonitor);

  // Set filename for the hdf5 ouput file
  std::string filename;
  if (ahandler.eventmode) {
    std::stringstream namestream;
    namestream << ahandler.outputfiletag << "_" << std::setw(4)
               << std::setfill('0') << ahandler.current_event << ".h5";
    filename = namestream.str();
  } else {
    filename = ahandler.outputfiletag + ".h5";
  }
  // Set file access which is append if we are restarting
  PetscFileMode file_access = FILE_MODE_WRITE;
  if (ahandler.restart) {
    file_access = FILE_MODE_APPEND;
  }
  // Open the file and create the measurement object
  Measurer measurer(model);

  // Creat the measurer output object, which will write the measurements
  // to the h5 file.
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

    // on separate frequency: measure coarsen solution
    if (ahandler.saveFrequencyCoarsen > 0 and
        steps % ahandler.saveFrequencyCoarsen == 0) {
      PetscLogEventBegin(measurements, 0, 0, 0, 0);
      measurer.measure_coarsen(&model->solution);
      if (rank == 0) {
        measurer_output->save_coarsen();
      }
      PetscLogEventEnd(measurements, 0, 0, 0, 0);
    }

    // on separate frequency: measure the topological charge Fourier readout
    if (ahandler.saveFrequencyTopcharge > 0 and
        steps % ahandler.saveFrequencyTopcharge == 0) {
      PetscLogEventBegin(measurements, 0, 0, 0, 0);
      measurer.measure_topcharge(&model->solution);
      if (rank == 0) {
        measurer_output->save_topcharge();
      }
      PetscLogEventEnd(measurements, 0, 0, 0, 0);
    }

    // Write the solution to tape if writeFrequency > 0. This is used for
    // plotting of the solution. It is normally not analyzed, or written.
    if (ahandler.writeFrequency > 0 and steps % ahandler.writeFrequency == 0) {
      PetscLogEventBegin(saving, 0, 0, 0, 0);
      std::ostringstream tString;
      // Set the precision to 2 digits after the decimal point filled with 0
      tString << std::fixed << std::setprecision(2) << "_t_" << atime.t();
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

  // Construct the stepper
  std::unique_ptr<Stepper> step;
  auto &etype = inputdata.ahandler.evolverType;
  if (etype == "PV2HBSplit23") {
    // Default is to include all steps
    PetscPrintf(
        PETSC_COMM_WORLD,
        "Using the default stepper PV2HBSplit23 with steps ABBABBABBC\n");
    step = std::make_unique<PV2HBSplit>(model, "ABBABBABBC", true, true, true);
  } else if (etype == "PV2HBSplitGeneral") {
    PetscPrintf(
        PETSC_COMM_WORLD,
        "Using the general stepper PV2HBSplit with steps from input file\n");
    nlohmann::json general_stepper = inputs["pv2hb_split_general"];
    std::string steps = general_stepper.value("steps", "ABBABBABBC");
    const bool ideal = general_stepper.value("include_ideal", true);
    const bool heatbath = general_stepper.value("include_heatbath", true);
    const bool diffuse = general_stepper.value("include_diffuse", true);
    step = std::make_unique<PV2HBSplit>(model, steps, ideal, heatbath, diffuse);
  } else if (etype == "ModelGDiffusionStep") {
    // Solves the diffusion equation
    PetscPrintf(PETSC_COMM_WORLD, "Using the ModelGDiffusionStep stepper with "
                                  "use_implicit step option\n");
    bool use_implicit_step = inputs["ModelGDiffusionStep"]["use_implicit_step"];
    if (use_implicit_step) {
      PetscPrintf(PETSC_COMM_WORLD, "Using implicit step\n");
      step = std::make_unique<ModelGDiffusionStep>(model, use_implicit_step);
    } else {
      PetscPrintf(PETSC_COMM_WORLD, "Using explicit step\n");
      step = std::make_unique<ModelGExplicitDiffusionStep>(model);
    }
  } else if (etype == "SuperSplitStep") {
    PetscPrintf(
        PETSC_COMM_WORLD,
        "Using the SuperSplitStep stepper with use_implicit_step option\n");
    bool use_implicit_step =
        inputs["SuperSplitStep"].value<bool>("use_implicit_step", false);
    std::string step_sequence =
        inputs["SuperSplitStep"].value<std::string>("step_sequence", "AB");
    PetscPrintf(PETSC_COMM_WORLD,
                "Using the SuperSplitStep stepper step_sequence %s \n",
                step_sequence.c_str());
    step = std::make_unique<SuperSplitStep>(model, step_sequence,
                                            use_implicit_step);
    ;
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "Unrecognized stepper type %s. Aborting...\n",
                etype.c_str());
    return;
  }

  auto &ahandler = model.data.ahandler;
  auto &atime = model.data.atime;
  int nevents = std::max(ahandler.nevents, 1);
  for (int i = 0; i < nevents; i++) {
    atime.reset();
    initialize_event(i, &model, inputs);
    run_event(i, &model, step.get(), inputs);
    ahandler.current_event++;
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
  nlohmann::json inputs;
  std::ifstream ifs(filename);
  if (ifs) {
    ifs >> inputs;
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "Unable to open input file %s. Aborting...\n",
                filename);
    return PetscFinalize();
  }
  PetscPrintf(PETSC_COMM_WORLD, "Current version: %s\n", gitversion);

  PetscPrintf(PETSC_COMM_WORLD, "Running SuperPions with input file %s\n",
              filename);

  int rank = 0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "Input parameters:\n" << inputs.dump(2) << std::endl;
  }

  Run(inputs);

  return PetscFinalize();
}
