#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <memory>
#include <vector>
#include <iomanip>

#include "ModelA.h"
#include "NoiseGenerator.h"
#include "Stepper.h"
#include "make_unique.h"
#include "plotter.h"

// Homemade parser
#include "parameterparser/parameterparser.h"

// Measurer, where the Petsc are included
#include "measurer.h"

struct ideal_data {
  double phi = 1;
  double NX;
  double deltaphi= 0.00001;
};

double ideal_fcn_flat(const double &x, const double &y, const double &z,
                           const int &L, void *params) {
  if (L == 3) {
	ideal_data *data = (ideal_data *)params;
      return data->phi*(1+data->deltaphi*
      exp(-0.1*(x-data->NX/2.)*(x-data->NX/2.)));
      //(1-cos(x*2*M_PI/data->NX));
      //(sin(x*2*M_PI/data->NX));
  	}
  else if (L == 0)  {
       ideal_data *data = (ideal_data *)params;
       return data->phi ; //some test data. This should not be modified by the evolutions
  	}
  else {
    	return 0;
       }
}

int main(int argc, char **argv) {
  // Initialization
  PetscErrorCode ierr;
  std::string help = "Model A. Call ./ModelA-Beuler.exe input=input.in with "
                     "input.in your input file.";

  ierr = PetscInitialize(&argc, &argv, (char *)0, help.c_str());
  if (ierr) {
    return ierr;
  }

  // Read the data form command line
  FCN::ParameterParser params(argc, argv);
  ModelAData inputdata(params);

  std::ofstream energy(inputdata.outputfiletag + "energy.txt");

  // If -quit flag, then just stop before we do anything but gather inputs.
  PetscBool quit = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-quit", &quit, NULL);
  CHKERRQ(ierr);
  if (quit) {
    return PetscFinalize();
  }

  // allocate the grid and initialize
  ModelA model(inputdata);
  ideal_data ideal;
  ideal.NX= (double) model.data.NX;
  model.initialize(ideal_fcn_flat,&ideal);

  // initialize the measurments and measure the initial condition
  Measurer measurer(&model);
  measurer.measure(&model.solution);

  // Initialize the stepper
  std::unique_ptr<Stepper> step;


  switch (model.data.evolverType) {
  case 1:
    step = make_unique<IdealLF>(model);
    break;
  case 2:
    step = make_unique<IdealPV2>(model);
    break;
  }
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  plotter plot(inputdata.outputfiletag + "_phi");
  //plot.plot(model.solution, "phi_0");
  PetscInt steps = 0;
    plot.settime(model.solution, steps,"phi");
    plot.dump(model.solution);
  // Thsi is the loop for the time step
  for (PetscReal time = model.data.initialtime; time < model.data.finaltime;
       time += model.data.deltat) {
    // Copy the solution
    VecCopy(model.solution, model.previoussolution);

    step->step(model.data.deltat);
    measurer.computeEnergy(model.data.deltat, model.data.evolverType == 1);

    if(myRank == 0){
      auto enComp = measurer.getEnergy();

      energy << std::setprecision(16) << time << " " << enComp[0]<< " " << enComp[1]<< " " << enComp[2]<< " " << enComp[3] << std::endl;
    }

    steps++;
      plot.update();
      plot.dump(model.solution);
    //std::string index= ;
    //std::string blabla = "phi_" + std::to_string(steps) ;
    //plot.plot(model.solution, "phi_" + std::to_string(steps));
  }

  // Destroy everything
  plot.finalize();
  step->finalize();
  measurer.finalize();
  model.finalize();
  energy.close();
  return PetscFinalize();
}
