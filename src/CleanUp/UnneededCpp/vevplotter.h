#ifndef VEVPLOTTER_H
#define VEVPLOTTER_H

#include <petscviewer.h>
#include <petscvec.h>
#include <petscdm.h>

class ModelA;

// Prints out the z=0 slice of the grid at regular time intervals into
// the hdf5 file named  $nametag_slices.h5.
//
// The usage is the following 
//
// ... create ModelA * model 
//
// VevPlotter plotter(model, time_per_analysis) ;
// 
// while(1) {
//
//    plotter.analyze(model)
//    ... step
//
// }
class VevPlotter {

private:
  unsigned int totalsteps = 0;
  unsigned int nsteps_per_analysis;

  PetscViewer H5viewer;
  Vec solution_slice;
  VecScatter scatter;

public:
  VevPlotter(ModelA *const model, const double &time_per_analysis=10.);
  ~VevPlotter() ;
  void analyze(ModelA *const model);
};


PetscErrorCode DMDAGetZSlice(DM da, Vec solution, PetscInt gp,
                             Vec *solution_slice, VecScatter *scatter) ;

#endif
