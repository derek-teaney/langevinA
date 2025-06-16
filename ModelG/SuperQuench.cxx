#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <vector>

#include "ModelA.h"
#include "NoiseGenerator.h"
#include "Stepper.h"
#include "make_unique.h"

// Measurer, where the Petsc are included
#include "measurer.h"
#include "measurer_output.h"


// call subroutines of ModelA that assign initial values to all fields (charges and phis)
void initialize_event(ModelA *model) {
  
  // call ModelA subroutine that initializes random spin configurations
  model->initialize_random_spins();
  // and subroutine that initializes gaussian random charges with normalization such that total
  // charge is zero
  model->initialize_gaussian_charges();


  // for former 'wave initial conditions' do instead:
  // model->initialize_wave_spins();
}


// This is the main loop of the program. It is a simple loop that steps the
// solution forward in time until the final time is reached. The data is analyzed
// and saved every saveFrequency steps.
void run_event(ModelA* const model,Stepper* const step) 
{

  const auto &ahandler = model->data.ahandler ;
  auto &atime = model->data.atime ;
  atime.reset() ;

  initialize_event(model) ;
  // Write the grid to a file so we can see the initial conditions
  // The "outputfiletag" is the base name of the file. The grid will be
  // written to a file with the name outputfiletag_save.h5. The _save.h5
  // is added by the write function.
  model->write(ahandler.outputfiletag + "_grid") ;

  // Set up logging for PETSc so we can find out how much time 
  // each part takes
  PetscInt steps = 0;
  PetscLogEvent measurements, stepmonitor, saving;
  PetscLogEventRegister("Measurements", 0, &measurements);
  PetscLogEventRegister("Saving the fields", 0, &saving);
  PetscLogEventRegister("Steps", 0, &stepmonitor);

  // Set filename for the hdf5 output of measurements
  std::string filename;
  if (ahandler.eventmode) {
    std::stringstream namestream;
    namestream << ahandler.outputfiletag << "_" << std::setw(4)
               << std::setfill('0') << ahandler.current_event << ".h5";
    filename = namestream.str();
  } else {
    filename = ahandler.outputfiletag + ".h5";
  }
  // Set file access mode for the hdf5 output of measurements
  PetscFileMode file_access = FILE_MODE_WRITE ;

  // Open the file and create the measurement object
  Measurer measurer(model) ;
  int rank = -1;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  std::unique_ptr<measurer_output_fasthdf5> measurer_output ;
  if (rank ==0) {
    measurer_output = std::make_unique<measurer_output_fasthdf5>(&measurer, filename, file_access);
  }

  // Start the loop
  const double tiny = 1.e-10;
  while (atime.t() < atime.tfinal() - tiny) {

    // measure the solution every saveFrequency
    if (steps % ahandler.saveFrequency == 0) {
      PetscLogEventBegin(measurements, 0, 0, 0, 0);
      measurer.measure(&model->solution);
      if (rank == 0) {
        measurer_output->save() ;
      }
      PetscPrintf(PETSC_COMM_WORLD,
                  "Event/Timestep %d/%d: step size = %g, time = %g, final = %g\n", ahandler.current_event, steps,
                  (double)atime.dt(), (double)atime.t(), (double)atime.tfinal());
      PetscLogEventEnd(measurements, 0, 0, 0, 0);
    }

    // Write the solution to tape if writeFrequency > 0. This is used for
    // plotting of the solution. It is normally not analyzed, or written.
    if(ahandler.writeFrequency > 0 and steps % ahandler.writeFrequency == 0) {
      PetscLogEventBegin(saving, 0, 0, 0, 0);
      std::ostringstream tString;
      tString << std::setprecision(4) <<"_t_" << atime.t();
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


int main(int argc, char **argv) {

  // Initialization of PETSc universe
  PetscErrorCode ierr;
  std::string help =
      "Usage: \n\n\t ./SuperQuench.exe -input input.json [options]\n\n";
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

  // Digest the inputs, some fields may be modified on ouptut.
  // If you need to add more parameters add them to ModelAData
  ModelAData inputdata(inputs);

  // If -quit flag, then just stop before we do anything but gather inputs.
  PetscBool quit = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-quit", &quit, NULL);
  CHKERRQ(ierr);
  if (quit) {
    return PetscFinalize();
  }

  // allocate the grid 
  ModelA model(inputdata);

  // Construct the stepper 
  std::unique_ptr<Stepper> step; 
  step = make_unique<IdealPV2>(model);

  auto &ahandler = model.data.ahandler ;
  if (ahandler.eventmode) { 
    for (int i = 0 ;  i < ahandler.nevents ; i++ )  {
      run_event(&model, step.get()) ;
      ahandler.current_event++ ;
    }
  } else  {
    run_event(&model, step.get()) ;
  }

  // Destroy everything
  step->finalize();
  model.finalize();
  return PetscFinalize();
}
