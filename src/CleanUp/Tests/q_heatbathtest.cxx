#include "ModelA.h"
#include "Stepper.h"
#include "make_unique.h"
#include "plotter.h"
// Homemade parser
#include "parameterparser/parameterparser.h"

void t1(const ModelAData &data, const double &nA0, const double &nB0,
        const std::string &filename);

struct diffusiontest_data {
  double w0 = 6.;
  double t = 0.;
  double A = 100000.;
  ModelAData &mdata;

  diffusiontest_data(ModelAData &in) : mdata(in) { ; }
};

double diffusiontest_fcn(const double &x, const double &y, const double &z,
                         const int &L, void *params) {
  if (L < ModelAData::Nphi) {
    return 0.;
  } else {
    diffusiontest_data *data = (diffusiontest_data *)params;

    ModelAData &mdata = data->mdata;

    double D = mdata.sigma() / mdata.chi;
    double w2 = pow(data->w0, 2) + 2. * D * data->t;

    double s = 0.;
    double N = mdata.NX / 2;
    static const int c = 2; // construct a periodic solution using c images
    for (int i = -c; i <= c; i++) {
      for (int j = -c; j <= c; j++) {
        for (int k = -c; k <= c; k++) {
          s += data->A *
               exp(-(pow(x - (2 * i + 1) * N, 2) + pow(y - (2 * j + 1) * N, 2) +
                     pow(z - (2 * k + 1) * N, 2)) /
                   (2 * w2)) /
               pow(2.0 * M_PI * w2, 3. / 2.);
        }
      }
    }
    return s;
  }
}

int main(int argc, char **argv) {
  // Initialization
  PetscErrorCode ierr;
  std::string help = "Tests a simple matrix inversion";
  ierr = PetscInitialize(&argc, &argv, (char *)0, help.c_str());
  if (ierr) {
    return ierr;
  }

  // Read the data form command line
  FCN::ParameterParser params(argc, argv);
  ModelAData inputdata(params);

  // allocate the grid and initialize
  ModelA model(inputdata);
  diffusiontest_data gauss(model.data);
  model.initialize(diffusiontest_fcn, &gauss);

  // Initialize the stepper
  // std::unique_ptr<Stepper> step = make_unique<ModelGChargeHB>(model);
  std::unique_ptr<Stepper> step = make_unique<ModelGChargeCN>(model);

  plotter plot(inputdata.outputfiletag + "_phi");
  plot.plot(model.solution, "phi_ic");

  plot.plotfcn(model.domain, model.solution, "phi_ic_sol", diffusiontest_fcn,
               &gauss);

  // Check the thermalization step
  if (model.rank == 0) {
    t1(inputdata, 0, 0, "t1_zero.out");
    t1(inputdata, 2, 1, "z1_21.out");

    // Check the listed faces
    for (int ixyz = 0; ixyz < 3; ixyz++) {
      for (int ieo = 0; ieo < 2; ieo++) {
        printf("%d %d %d %d \n", g_face_cases[ixyz][ieo].eoA,
               g_face_cases[ixyz][ieo].iB, g_face_cases[ixyz][ieo].jB,
               g_face_cases[ixyz][ieo].kB);
      }
    }
  }

  PetscInt steps = 1;
  // Thsi is the loop for the time step
  PetscReal time;
  for (time = model.data.initialtime; time < model.data.finaltime;
       time += model.data.deltat) {
    VecCopy(model.solution, model.previoussolution);
    step->step(model.data.deltat);

    if (steps % model.data.saveFrequency == 0) {
      // Print some information to not get bored during the running:
      PetscPrintf(PETSC_COMM_WORLD, "Timestep %D: step size = %g, time = %g\n",
                  steps, (double)model.data.deltat, (double)time);
    }
    steps++;
  }

  plot.plot(model.solution, "phi_final");
  gauss.t = time;
  plot.plotfcn(model.domain, model.solution, "phi_final_sol", diffusiontest_fcn,
               &gauss);

  step->finalize();
  plot.finalize();
  model.finalize();

  return PetscFinalize();
}

void t1(const ModelAData &data, const double &nA0, const double &nB0,
        const std::string &filename) {
  o4_stepper_monitor monitor;
  FILE *fp = fopen(filename.c_str(), "w");

  const double rms = sqrt(2. * data.deltat * data.sigma());
  const int N = 500000;
  double nA = nA0;
  double nB = nB0;
  for (int i = 0; i < N; i++) {
    double q = modelg_update_charge_pair(data.chi, rms, nA, nB, monitor);
    nA -= q;
    nB += q;
    fprintf(fp, "%16.8e %16.8e %15.8e \n", nA, nB, q);
  }
  fclose(fp);
}
