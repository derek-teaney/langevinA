#include "vevplotter.h"
#include "ModelA.h"
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <petscviewerhdf5.h>

PetscErrorCode DMDAGetZSlice(DM da, Vec solution, PetscInt gp,
                             Vec *solution_slice, VecScatter *scatter) {

  // Get info about the dimensions of the grid
  PetscInt dim, M, N, P, dof, s;
  DMDAGetInfo(da, &dim, &M, &N, &P, NULL, NULL, NULL, &dof, &s, NULL, NULL,
              NULL, NULL);

  // This code is only designed for 3D domains
  if (dim != 3) {
    SETERRQ(PetscObjectComm((PetscObject)da), PETSC_ERR_SUP,
            "Cannot get slice from 1d or 2d DMDA");
  }

  // Get the rank. The new vector is created on rank 0
  PetscMPIInt rank;
  int ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)da), &rank);
  CHKERRMPI(ierr);

  // AO object is used to convert natural ordering of indices of grid  to PETSc
  // ordering. The AO object is constructed by the DMDA grid object
  AO ao;
  DMDAGetAO(da, &ao);

  IS is;
  if (rank == 0) {
    // Size of the new 2d array
    PetscInt dofMN = dof * M * N;

    // Create the indices array holding the natural (Application Ordered)
    // indices for the slice.
    PetscInt *indices;
    ierr = PetscMalloc1(dofMN, &indices); // indices is owned by IS object below
    CHKERRQ(ierr);

    // indices[0] is the starting index of the z slize, the subsequenty
    // elements are sequentual
    indices[0] = gp * dofMN;
    for (int i = 1; i < dofMN; i++) {
      indices[i] = indices[i - 1] + 1;
    }
    ierr = AOApplicationToPetsc(ao, dofMN, indices);
    CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_SELF, dofMN, indices, PETSC_OWN_POINTER,
                           &is);
    CHKERRQ(ierr);

    // Create  a 2d domain for the slice and the corresponding vector
    DM da2d;
    DMDACreate2d(MPI_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_STAR, M, N, PETSC_DECIDE, PETSC_DECIDE, dof, s,
                 NULL, NULL, &da2d);
    DMSetUp(da2d);
    DMCreateGlobalVector(da2d, solution_slice);
    DMDestroy(&da2d);

  } else {
    ierr = VecCreateSeq(PETSC_COMM_SELF, 0, solution_slice);
    CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_SELF, 0, NULL, PETSC_COPY_VALUES, &is);
    CHKERRQ(ierr);
  }

  // Vec is any vector with the shape of the vector we scatter from
  // While solution_slice is the shape we will scatter to.  is is the indices
  // we scatter from.
  ierr = VecScatterCreate(solution, is, *solution_slice, NULL, scatter);
  CHKERRQ(ierr);
  ierr = ISDestroy(&is);
  CHKERRQ(ierr);

  return (0);
}

/////////////////////////////////////////////////////////////////////////
VevPlotter::VevPlotter(ModelA *const model, const double &time_per_analysis)
    : totalsteps(0), nsteps_per_analysis(1) {

  nsteps_per_analysis = static_cast<unsigned int>(time_per_analysis/model->data.atime.dt()) ;

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Initializaing Slice %d\n", rank)  ;
  // Set the datatransfer mechanism
  PetscInt zslice = 0;
  DMDAGetZSlice(model->domain, model->solution, 0, &solution_slice, &scatter);

  // Set the name of the dataset in the hdf5 file
  PetscInt ierr =
      PetscObjectSetName((PetscObject)solution_slice, "solution_slices");
  CHKERRV(ierr);

  // Open the hdf5 file for writing
  if (rank == 0) {
    std::string name = model->data.ahandler.outputfiletag + "_slices.h5";
    ierr = PetscViewerHDF5Open(PETSC_COMM_SELF, name.c_str(), FILE_MODE_WRITE,
                               &H5viewer);
    CHKERRV(ierr);

    ierr = PetscViewerHDF5SetTimestep(H5viewer, 0);
    CHKERRV(ierr);
  }
}

void VevPlotter::analyze(ModelA *const model) {

  // It is not time to do the analysis, just increment the counter
  if (totalsteps % nsteps_per_analysis != 0) {
    totalsteps++;
    return;
  }
  // It is time to do the analysis increment the step counter and do the
  // analysis
  totalsteps++;

  // Do the actual transfer of data
  VecScatterBegin(scatter, model->solution, solution_slice, INSERT_VALUES,
                  SCATTER_FORWARD);
  VecScatterEnd(scatter, model->solution, solution_slice, INSERT_VALUES,
                SCATTER_FORWARD);

  // Write the slice to the hdf5 file
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    VecView(solution_slice, H5viewer);
    PetscViewerHDF5IncrementTimestep(H5viewer);
    return;
  } else {
    return;
  }
}

VevPlotter::~VevPlotter() {

  VecDestroy(&solution_slice);
  VecScatterDestroy(&scatter);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    PetscPrintf(PETSC_COMM_WORLD, "Tearing down HDF5 %d\n", rank)  ;
    PetscViewerDestroy(&H5viewer);
    PetscPrintf(PETSC_COMM_WORLD, "Completed Tear down\n") ;
  }
}
