#ifndef HDF5_NTUPLE_H
#define HDF5_NTUPLE_H

#include "hdf5.h"
#include <array>
#include <vector>
#include <string>


/*  This gives the basic usage 
int main(int argc, char **argv) {

  hid_t file_id;
  file_id = H5Fcreate("work.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  static const hsize_t N = 8;
  std::array<double, N> datain = {1, 2, 3, 4, 5, 6, 7, 8};

  // Create a 2x4  2D array data set on the file which will be
  // determined at each time step. So on time it is three dimensional
  //
  // (unlimitted, 2, 4)
  //
  std::array<size_t, 2> NN2 = {2,4} ;
  ntuple<2> nt2(NN2, "phi2", file_id);

  for (size_t j = 0; j < 3; j++) {

    for (int i0=0; i0 < NN2[0] ; i0++) {
    for (int i1=0; i1 < NN2[1] ; i1++) {
      size_t k = nt2.at({i0, i1}) ; //Mapping the 2D array to a single index
      nt2.row[k] = datain[k] ;  // fill up the buffer row
    }
    }

    nt2.fill(); // write the buffer to the file

  }

  // Close the ntuple
  nt2.close();

  // close the filespace
  H5Fclose(file_id);
  return 0;
}
*/


template<std::size_t rank>
class ntuple {

public:
  std::vector<double> row ;

  size_t at(const std::array<int, rank> &idex) {
    size_t offset = 0;
    for(size_t i=0; i < rank; i++) offset += sizes[i]*idex[i] ;
    return offset ;
  }

  ntuple(const std::array<size_t, rank> &n, const std::string &nm,
         const hid_t &file_id) : closed(false) {

    size_t count =1 ;
    for (size_t i = 0; i < rank; i++) {
      N[i] = n[i];
      count *= n[i] ;
    }
    row.resize(count) ;

    sizes[rank-1] = 1 ;
    // gotcha: use int not size_t can go negative
    for (int i = static_cast<int>(rank)-2; i >=0; i--) {
       sizes[i] = N[i+1]*sizes[i+1] ;
    }
    
    // Open the dataset
    htri_t exists = H5Lexists(file_id, nm.c_str(), H5L_TYPE_HARD) ;
    if (!exists) {
      create(file_id, nm) ;
    } else {
      open(file_id, nm) ;
    }

    // Memory of space
    hsize_t line_dims[rank+1];
    pack(line_dims, 1, N);
    memory_space = H5Screate_simple(rank+1, line_dims, NULL);
  }

  ~ntuple() {
    if (!closed) close();
  }

  void close() {
    if (closed) return ;

    // Close the filespace
    H5Sclose(memory_space);

    // and the dataset
    H5Dclose(dataset);

    // and the dataspace
    H5Sclose(dataspace);

    closed = true ;
  }

  void fill() {
    // Extend by a row
    dims[0]++;
    herr_t status = H5Dextend(dataset, dims);

    // Get the space in the file
    hid_t filespace = H5Dget_space(dataset);

    // set the offset and count
    hsize_t offset[rank+1] = {};
    offset[0] = dims[0] - 1;
    hsize_t count[rank+1];
    pack(count, 1, N);

    // Get the hyperslab
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count,
                                 NULL);

    // Write the data
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memory_space, filespace,
                      H5P_DEFAULT, row.data());

    // Close the filespace
    H5Sclose(filespace);
  }

private:

  // number of columns
  hsize_t N[rank];

  // memory
  bool closed;
  hid_t dataspace;
  hid_t dataset;
  hid_t memory_space;

  // dimensions of the stored array
  hsize_t dims[rank+1];

  // helper
  void pack(hsize_t d[], hsize_t d0, hsize_t n[]) {
    d[0] = d0 ;
    for (size_t i = 0; i< rank; i++) {
      d[i+1] = n[i] ;
    }
  }

  size_t sizes[rank] ;

  void create(hid_t file_id, const std::string &nm)  {
    // Set the  array dimensions
    pack(dims, 0, N);
    hsize_t maxdims[rank+1];
    pack(maxdims, H5S_UNLIMITED, N);

    // Create the dataspace
    dataspace = H5Screate_simple(rank+1, dims, maxdims);

    // set chunking properties of the dataset
    hid_t chunk_params = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dims[rank+1];
    pack(chunk_dims, 2, N);
    herr_t status = H5Pset_chunk(chunk_params, rank+1, chunk_dims);

    // Open the dataset
    dataset = H5Dcreate(file_id, nm.c_str(), H5T_NATIVE_DOUBLE, dataspace,
                        H5P_DEFAULT, chunk_params, H5P_DEFAULT);

    H5Pclose(chunk_params);
  }

  void open(hid_t file_id, const std::string &nm) {

    dataset = H5Dopen(file_id,  nm.c_str(), H5P_DEFAULT) ;
    dataspace = H5Dget_space(dataset);    

    // This line reads the dimensions of the array from tape
    hsize_t maxdims[rank+1] ;
    herr_t status  = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

    // check that the dimensions agree
    for (size_t i=0; i < rank; i++) {
      if (N[i] != dims[i+1]) {
        throw ("ntuple::ntuple Error dataset dimensions do not conform") ;
      }
    }
  }
  
};

#endif
