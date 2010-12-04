// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      disk_FileHdf5.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:11:36 PST 2008
/// @brief     Implementation of the FileHdf5 class

#include <string>

#include "error.hpp"
#include "disk.hpp"
 
//----------------------------------------------------------------------

FileHdf5::FileHdf5()
/**
 */
 :file_(0),
  file_name_(),
  file_mode_(),
  is_file_open_(false),
  dataset_(0),
  dataset_name_(),
  dataspace_(0)
{
#ifdef CONFIG_USE_HDF5
#endif
}

//----------------------------------------------------------------------
    
int FileHdf5::file_open  (std::string name, std::string mode)
/**
 * @param name  Name of the file to create or open
 * @param mode  How the file is to be created or opened:
 *              - "r": read-only: file must exist.
 *              - "w": write-only: any existing file will be truncated.
 * @return      True iff opening the file was successful
 */
{
#ifdef CONFIG_USE_HDF5
  name = name + ".hdf5";

  if (is_file_open_) {

    char warning_message[ERROR_MESSAGE_LENGTH];
    sprintf (warning_message,"Attempting to open an open file %s",name.c_str());
    WARNING_MESSAGE("FileHdf5::file_open",warning_message);

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
      ERROR_MESSAGE("FileHdf5::file_open",error_message);
    }

    if (file_ >= 0) {
      is_file_open_ = true;
    } else {
      char warning_message[ERROR_MESSAGE_LENGTH];
      sprintf (warning_message,
	       "Return value %d opening file %s",
	       file_,file_name_.c_str());
      WARNING_MESSAGE("FileHdf5::file_open",warning_message);
    }
  }

  return is_file_open_;
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

void FileHdf5::file_close ()
/**
 */
{
#ifdef CONFIG_USE_HDF5
  if (! is_file_open_) {
    char warning_message[ERROR_MESSAGE_LENGTH];
    sprintf (warning_message,
	     "Attempting to close a closed file %s",
	     this->file_name_.c_str());
    WARNING_MESSAGE("FileHdf5::file_close",warning_message);
  } else {
    int retval = H5Fclose (file_);
    if (retval >= 0) {
      is_file_open_ = false;
    } else {
      char warning_message[ERROR_MESSAGE_LENGTH];
      sprintf (warning_message,
	       "Return value %d closing file %s",
	       retval,file_name_.c_str());
      WARNING_MESSAGE("FileHdf5::file_close",warning_message);
    }
  }
#endif
}

//----------------------------------------------------------------------

void FileHdf5::group_open (std::string name)
/**
 */
{
#ifdef CONFIG_USE_HDF5
#endif
}

//----------------------------------------------------------------------

void FileHdf5::group_close ()
/**
 */
{
#ifdef CONFIG_USE_HDF5
#endif
}

//----------------------------------------------------------------------

void FileHdf5::dataset_open_write (std::string name, 
			       int nx,  int ny,  int nz)
			       //			       int ix0, int iy0, int iz0,  /// unused
			       //			       int mx,  int my,  int mz)   /// unused
{
#ifdef CONFIG_USE_HDF5
  if (file_mode_ != "w") {

    char error_message[ERROR_MESSAGE_LENGTH];
    sprintf (error_message, "Expecting dataset_open_read()");
    ERROR_MESSAGE("FileHdf5::dataset_open_write",error_message);

  } else {

    // Open the dataset for writing

    // determine the dimension d

    int d = 3;
    
    if (nz == 1) { d--;  if (ny == 1) { d--; }}

    ASSERT ("dataset_open_write","d is out of range",1 <= d && d <= 3);

    hsize_t n[3];

    // Transpose Enzo -> HDF5 row/column ordering

//     if (d == 1) {
//       n[0] = nx;
//     } else if (d == 2) {
//       n[0] = ny; // swap x, y
//       n[1] = nx;
//     } else if (d == 3) {
//       n[0] = nz; // swap x, y, z
//       n[1] = ny;
//       n[2] = nx;
//     }
      
    n[0] = nx;
    n[1] = ny;
    n[2] = nz;

    dataspace_ = H5Screate_simple (d,n,0);
    
    dataset_   = H5Dcreate( file_, 
			    name.c_str(), 
			    datatype, 
			    dataspace_,  
			    H5P_DEFAULT );

    if (dataset_ >= 0) {

      dataset_name_ = name;
      
    } else {

      char warning_message[ERROR_MESSAGE_LENGTH];
      sprintf (warning_message,
	       "Return value %d opening dataset %s",
	       dataset_,name.c_str());
      WARNING_MESSAGE("FileHdf5::dataset_open_write",warning_message);
    }
  }
#endif
}

//----------------------------------------------------------------------

void FileHdf5::dataset_open_read (std::string name, 
			      int *       nx,
			      int *       ny,
			      int *       nz)

/**
 */
{
#ifdef CONFIG_USE_HDF5
  if (file_mode_ != "r") {

    char error_message[ERROR_MESSAGE_LENGTH];
    sprintf (error_message, "Expecting dataset_open_write()");
    ERROR_MESSAGE("FileHdf5::dataset_open_read",error_message);

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
	ERROR_MESSAGE("FileHdf5::dataset_open",error_message);
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
      WARNING_MESSAGE("FileHdf5::dataset_open_read",warning_message);
    }
  }
#endif
}

//----------------------------------------------------------------------

// void FileHdf5::dataset_set_hyperslab ()
// {
//   hsize_t offset[3] = {ix0,iy0,iz0};
//   hsize_t count[3]  = {mx,my,mz};
//   status_ = H5Sselect_hyperslab 
//     (dataspace_, H5S_SELECT_SET,  offset, NULL, count, NULL);
// }

//----------------------------------------------------------------------

void FileHdf5::dataset_close ()
/**
 */
{
#ifdef CONFIG_USE_HDF5
  H5Dclose( dataset_);
#endif
}

//----------------------------------------------------------------------

void FileHdf5::read  (Scalar * buffer)
/**
 */
{
#ifdef CONFIG_USE_HDF5
  H5Dread (dataset_, datatype, dataspace_, H5S_ALL, H5P_DEFAULT, buffer);
#endif
}

//----------------------------------------------------------------------

void FileHdf5::write (Scalar * buffer,)
/**
 */
{
#ifdef CONFIG_USE_HDF5
  H5Dwrite (dataset_, datatype, dataspace_, H5S_ALL, H5P_DEFAULT, buffer);
#endif
}


  //----------------------------------------------------------------------

#ifdef CONFIG_USE_HDF5
  int precision_hdf5_(enum precision_enum precision)
  {
#   define SCALAR_HDF5    H5T_NATIVE_FLOAT
#   define SCALAR_HDF5    H5T_NATIVE_DOUBLE
    // NATIVE   X    FLOAT DOUBLE LDOUBLE
    // IEEE     X    F32BE F64BE
    // STD     B16BE B32BE B64BE 
    switch (precision) {
    case precision_unknown:
      return 0;
      break;
    case precision_default:
      return precision_hdf5(default_precision);
      break;
    case precision_half:
      return H5T_STD_B16BE; // 16-bit bit field (big-endian)
      break;
    case precision_single:
      return H5T_IEEE_F32BE; // 32-bit IEEE float (big-endian)
      break;
    case precision_double:
      return H5T_IEEE_F64BE; // 64-bit IEEE float (big-endian)
      break;
    case precision_extended80:
      return (sizeof(long double)==10);
      break;
    case precision_extended96:
      return (sizeof(long double)==12);
      break;
    case precision_quadruple:
      return (sizeof(long double)==16);
      break;
    default:
      return 0;
      break;
    }
  }
