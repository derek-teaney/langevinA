#include <memory>

#include "ModelA.h"
#include "NoiseGenerator.h"
#include "Stepper.h"
#include "plotter.h"

// Homemade parser
#include "parameterparser/parameterparser.h"

// Measurer, where the Petsc are included
#include "measurer.h"

struct smoothtest_data {
  double N = 32;
  double sigma = 4.;
};

double smoothtest_fcn_flat(const double &x, const double &y, const double &z,
                           const int &L, void *params) {
  if (L < ModelAData::Nphi) {
    smoothtest_data *data = (smoothtest_data *)params;
    double sigma2 = pow(data->sigma + L * 0.2, 2);
    double N = data->N / 2;
    return exp(-(pow(x - N + 0.5, 2) + pow(y - N - 1.0, 2) +
                 pow(z - N + 2, 2)) /
               (2 * sigma2)) /
           pow(2.0 * M_PI * sigma2, 3. / 2.);
  } else {
    return L; // some test data. This should not be modified by the evolutions
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

  // allocate the grid and initialize
  smoothtest_data gauss;
  ModelA model(inputdata);
  model.initialize(smoothtest_fcn_flat, &gauss);

  plotter plot(inputdata.outputfiletag + "_phi");
  plot.plot(model.solution, "phi_ic");

  // Initialize the stepper
  std::unique_ptr<Stepper> step;

  // Initialize the data steppers w/out noise
  switch (model.data.evolverType) {
  case 1:
    step = std::make_unique<BackwardEuler>(model, false);
    break;
  case 2:
    step = std::make_unique<ForwardEuler>(model, false);
    break;
  case 3:
    step = std::make_unique<SemiImplicitBEuler>(model, false);
    break;
  }

  PetscInt steps = 1;
  // Thsi is the loop for the time step
  for (PetscReal time = model.data.initialtime; time < model.data.finaltime;
       time += model.data.deltat) {
    VecCopy(model.solution, model.previoussolution);
    step->step(model.data.deltat);
    steps++;
  }

  plot.plot(model.solution, "phi_final");
  // Destroy everything
  plot.finalize();
  step->finalize();
  model.finalize();
  return PetscFinalize();
}
