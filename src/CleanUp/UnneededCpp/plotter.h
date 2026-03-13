#ifndef MODELA_PLOTTER
#define MODELA_PLOTTER

// This simple class is used to printout the fields phi at regular intervals
// into an hdf5 file
class plotter {

public:
  plotter(const std::string &filename) {
    std::string name = filename + ".h5";
    PetscViewerHDF5Open(PETSC_COMM_WORLD, name.c_str(), FILE_MODE_WRITE,
                        &H5viewer);
  }

    void settime(Vec V,int time,const std::string &name){
        PetscObjectSetName((PetscObject)V, name.c_str());
        PetscViewerHDF5SetTimestep(H5viewer, time);
        PetscViewerSetFromOptions(H5viewer);
    }
    
    void update() {
      PetscViewerHDF5IncrementTimestep(H5viewer);
    }
    
    void dump(Vec V) {
      VecView(V, H5viewer);
    }
  
    
  void plot(Vec V, const std::string &name) {
    PetscObjectSetName((PetscObject)V, name.c_str());
    VecView(V, H5viewer);
  }

  void plotfcn(DM domain, Vec v, const std::string &name,
               double (*fcn)(const double &x, const double &y, const double &z,
                             const int &l, void *params),
               void *params) {
    Vec f;
    VecDuplicate(v, &f);

    PetscScalar ****u;
    DMDAVecGetArrayDOF(domain, f, &u);

    PetscInt i, j, k, l, xstart, ystart, zstart, xdimension, ydimension,
        zdimension;

    PetscInt dof;
    DMDAGetInfo(domain, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &dof, NULL,
                NULL, NULL, NULL, NULL);

    DMDAGetCorners(domain, &xstart, &ystart, &zstart, &xdimension, &ydimension,
                   &zdimension);

    for (k = zstart; k < zstart + zdimension; k++) {
      for (j = ystart; j < ystart + ydimension; j++) {
        for (i = xstart; i < xstart + xdimension; i++) {
          for (l = 0; l < dof; l++) {
            u[k][j][i][l] = fcn(k, j, i, l, params);
          }
        }
      }
    }
    DMDAVecRestoreArrayDOF(domain, f, &u);
    plot(f, name);
    VecDestroy(&f);
  }

  void finalize() { PetscViewerDestroy(&H5viewer); }

private:
  PetscViewer H5viewer;
};

#endif
