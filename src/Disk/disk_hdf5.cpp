//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      hdf5.cpp
 * @brief     Implementation of the Hdf5 class
 * @author    James Bordner
 * @date      Thu Feb 21 16:11:36 PST 2008
 *
 * $Id$
 *
 *********************************************************************
 */

#include <hdf5.h>

#include "error.hpp"

#include "disk_hdf5.hpp"
 
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
  dataspace_(0),
  datatype_(SCALAR_HDF5)
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

    char warning_message[ERROR_MESSAGE_LENGTH];
    sprintf (warning_message,"Attempting to open an open file %s",name.c_str());
    WARNING_MESSAGE("Hdf5::file_open",warning_message);

  } else {

    file_name_ = name;
    file_mode_ = mode;

    if (mode == "r") {
      file_ = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    } else if (mode == "w") {
      file_ = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    } else {
      char error_message[ERROR_MESSAGE_LENGTH];
      sprintf (error_message,"Unrecognized mode: %s",mode.c_str());
      ERROR_MESSAGE("Hdf5::file_open",error_message);
    }

    if (file_ >= 0) {
      is_file_open_ = true;
    } else {
      char warning_message[ERROR_MESSAGE_LENGTH];
      sprintf (warning_message,
	       "Return value %d opening file %s",
	       file_,file_name_.c_str());
      WARNING_MESSAGE("Hdf5::file_open",warning_message);
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
    char warning_message[ERROR_MESSAGE_LENGTH];
    sprintf (warning_message,
	     "Attempting to close a closed file %s",
	     this->file_name_.c_str());
    WARNING_MESSAGE("Hdf5::file_close",warning_message);
  } else {
    int retval = H5Fclose (file_);
    if (retval >= 0) {
      is_file_open_ = false;
    } else {
      char warning_message[ERROR_MESSAGE_LENGTH];
      sprintf (warning_message,
	       "Return value %d closing file %s",
	       retval,file_name_.c_str());
      WARNING_MESSAGE("Hdf5::file_close",warning_message);
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

void Hdf5::dataset_open_write (std::string name, 
			       int nx,
			       int ny,
			       int nz)
{
  if (file_mode_ != "w") {

    char error_message[ERROR_MESSAGE_LENGTH];
    sprintf (error_message, "Expecting dataset_open_read()");
    ERROR_MESSAGE("Hdf5::dataset_open_write",error_message);

  } else {

    // Open the dataset for writing

    // determine the dimension d

    int d = 3;
    
    if (nz == 1) { // <= 2D
      d--;
      if (ny == 1) { // <= 1D
	d--;
      }
    }

    assert (1 <= d && d <= 3);

    // Transpose Enzo -> HDF5 row/column ordering

    hsize_t n[3];

    if (d == 1) {
      n[0] = nx;
    } else if (d == 2) {
      n[0] = ny; // swap x, y
      n[1] = nx;
    } else if (d == 3) {
      n[0] = nz; // swap x, y, z
      n[1] = ny;
      n[2] = nx;
    }
      
    dataspace_ = H5Screate_simple (d,n,0);
    
    dataset_   = H5Dcreate( file_, 
			    name.c_str(), 
			    datatype_, 
			    dataspace_,  
			    H5P_DEFAULT );

    if (dataset_ >= 0) {

      dataset_name_ = name;
      
    } else {

      char warning_message[ERROR_MESSAGE_LENGTH];
      sprintf (warning_message,
	       "Return value %d opening dataset %s",
	       dataset_,name.c_str());
      WARNING_MESSAGE("Hdf5::dataset_open_write",warning_message);
    }
  }
}

//----------------------------------------------------------------------

void Hdf5::dataset_open_read (std::string name, 
			      int *       nx,
			      int *       ny,
			      int *       nz)

/**
 */
{
  if (file_mode_ != "r") {

    char error_message[ERROR_MESSAGE_LENGTH];
    sprintf (error_message, "Expecting dataset_open_write()");
    ERROR_MESSAGE("Hdf5::dataset_open_read",error_message);

  } else {

    int d;        // Dataset dimension
    hsize_t n[3]; // Dataset size

    if (file_mode_ == "r") {

      // Open the dataset for reading

      // Open the dataset
      dataset_ = H5Dopen( file_, name.c_str());
      // Get the size of dataset
      dataspace_ = H5Dget_space (dataset_);
      d = H5Sget_simple_extent_ndims(dataspace_);
      if (d > 3) {
	char error_message[ERROR_MESSAGE_LENGTH];
	sprintf (error_message, "Dataset has too many dimensions %d",d);
	ERROR_MESSAGE("Hdf5::dataset_open",error_message);
      }

      H5Sget_simple_extent_dims(dataspace_,n,0);

      // Set the array size accordingly
    
      *nx = n[0];
      *ny = (d > 1) ? n[1] : 1;
      *nz = (d > 2) ? n[2] : 1;

    }

    if (dataset_ >= 0) {

      dataset_name_ = name;
      
    } else {

      char warning_message[ERROR_MESSAGE_LENGTH];
      sprintf (warning_message,
	       "Return value %d opening dataset %s",
	       dataset_,name.c_str());
      WARNING_MESSAGE("Hdf5::dataset_open_read",warning_message);
    }
  }
}

//----------------------------------------------------------------------

void Hdf5::dataset_close ()
/**
 */
{
  H5Dclose( dataset_);
}

//----------------------------------------------------------------------

void Hdf5::read  (Scalar * buffer)
/**
 */
{
  H5Dread (dataset_, datatype_, dataspace_, H5S_ALL, H5P_DEFAULT, buffer);
}

//----------------------------------------------------------------------

void Hdf5::write (Scalar * buffer)
/**
 */
{
  H5Dwrite (dataset_, datatype_, dataspace_, H5S_ALL, H5P_DEFAULT, buffer);
}

