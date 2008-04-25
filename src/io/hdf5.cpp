/** 
 *********************************************************************
 *
 * @file      hdf5.cpp
 * @brief     Implementation of the Hdf5 class
 * @author    James Bordner
 * @date      Thu Feb 21 16:11:36 PST 2008
 *
 *********************************************************************
 */

#include <string>

#include "hdf5.h"

#include "error.hpp"
#include "scalar.hpp"
#include "array.hpp"
#include "hdf5.hpp"
 
//----------------------------------------------------------------------

Hdf5::Hdf5()
/**
 */
 :file_(0),
  file_name_(),
  file_mode_(),
  is_file_open_(false),
  dataset_(0),
  dataset_name_(),
  is_dataset_open_(false),
  dataspace_(0),
  datatype_(H5T_NATIVE_DOUBLE)
{
}

//----------------------------------------------------------------------
    
int Hdf5::file_open  (std::string name, std::string mode)
/**
 * @param name  Name of the file to create or open
 * @param mode  How the file is to be created or opened:
 *              - "r": read-only: file must exist.
 *              - "w": write-only: any existing file will be truncated.
 * @return      True iff opening the file was successful
 */
{
  if (is_file_open_) {

    sprintf (warning_message,"Attempting to open an open file %s",name.c_str());
    WARNING_MESSAGE("Hdf5::file_open");

  } else {

    file_name_ = name;
    file_mode_ = mode;

    if (mode == "r") {
      file_ = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    } else if (mode == "w") {
      file_ = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    } else {
      sprintf (error_message,"Unrecognized mode: %s",mode.c_str());
      ERROR_MESSAGE("Hdf5::file_open");
    }

    if (file_ >= 0) {
      is_file_open_ = true;
    } else {
      sprintf (warning_message,
	       "Return value %d opening file %s",
	       file_,file_name_.c_str());
      WARNING_MESSAGE("Hdf5::file_open");
    }
  }

  return is_file_open_;

}

//----------------------------------------------------------------------

void Hdf5::file_close ()
/**
 */
{
  if (! is_file_open_) {
    sprintf (warning_message,
	     "Attempting to close a closed file %s",
	     this->file_name_.c_str());
    WARNING_MESSAGE("Hdf5::file_close");
  } else {
    int retval = H5Fclose (file_);
    if (retval >= 0) {
      is_file_open_ = false;
    } else {
      sprintf (warning_message,
	       "Return value %d closing file %s",
	       retval,file_name_.c_str());
      WARNING_MESSAGE("Hdf5::file_close");
    }
  }
}

//----------------------------------------------------------------------

void Hdf5::group_open (std::string name)
/**
 */
{
}

//----------------------------------------------------------------------

void Hdf5::group_close ()
/**
 */
{
}

//----------------------------------------------------------------------

void Hdf5::dataset_open (std::string name, Array & array)
/**
 */
{
  int d;        // Dataset dimension
  hsize_t n[3]; // Dataset size

  if (file_mode_ == "r") {
    // Open the dataset
    dataset_ = H5Dopen( file_, name.c_str());
    // Get the size of dataset
    dataspace_ = H5Dget_space (dataset_);
    d = H5Sget_simple_extent_ndims(dataspace_);
    if (d > 3) {
      sprintf (error_message, "Dataset has too many dimensions %d\n",d);
      ERROR_MESSAGE("Hdf5::dataset_open");
    }
    H5Sget_simple_extent_dims(dataspace_,n,0);
    // Set the array size accordingly
    if (d == 1) array.resize(n[0]);
    if (d == 2) array.resize(n[0],n[1]);
    if (d == 3) array.resize(n[0],n[1],n[2]);
  } else if (file_mode_ == "w") {

    int m[3];
    array.size(&m[0],&m[1],&m[2]);
    d = 3;
    if (n[2] == 1) {
      d--;
      if (n[1] == 1) {
	d--;
      }
    }

    n[0] = m[0];
    n[1] = m[1];
    n[2] = m[2];
    dataspace_ = H5Screate_simple (d,n,0);
    
    dataset_   = H5Dcreate( file_, name.c_str(), datatype_, dataspace_,  H5P_DEFAULT);
  }

  if (dataset_ >= 0) {

    is_dataset_open_ = true;
    dataset_name_ = name;
      
  } else {
    sprintf (warning_message,
	     "Return value %d opening dataset %s",
	     dataset_,name.c_str());
    WARNING_MESSAGE("Hdf5::dataset_open");
  }
}

//----------------------------------------------------------------------

void Hdf5::dataset_close ()
/**
 */
{
}

//----------------------------------------------------------------------

void Hdf5::read  (Array & array)
/**
 */
{
  Scalar * buffer = array.values();

  H5Dread (dataset_, datatype_, dataspace_, H5S_ALL, H5P_DEFAULT, buffer);
}

//----------------------------------------------------------------------

void Hdf5::write (Array & array)
/**
 */
{
  Scalar * buffer = array.values();

  H5Dwrite (dataset_, datatype_, dataspace_, H5S_ALL, H5P_DEFAULT, buffer);
}

