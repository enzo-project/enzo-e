// See LICENSE_CELLO file for license and copyright information

/// @file      disk_FileHdf5.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:11:36 PST 2008
/// @brief     Implementation of the FileHdf5 class

#include "cello.hpp"

#include "disk.hpp"
 
//----------------------------------------------------------------------

FileHdf5::FileHdf5
(
 std::string path,
 std::string name,
 std::string mode
 ) throw()
  : File(path,name,mode),
    file_id_(0),
    data_set_id_(0),
    data_space_id_(0),
    status_id_(0),
    data_name_(""),
    data_type_(scalar_type_unknown),
    data_rank_(0),
    is_file_open_(false),
    is_data_open_(false)
{
  const int rank_max = 5;
  for (int i=0; i<rank_max; i++) {
    data_size_[i] = 0;
  }
}

//----------------------------------------------------------------------
    
void FileHdf5::open () throw()
/**
 * @param name  Name of the file to create or open
 * @param mode  How the file is to be created or opened:
 *              - "r": read-only: file must exist.
 *              - "w": write-only: any existing file will be truncated.
 * @return      True iff opening the file was successful
 */
{

  if (is_file_open_) {

    char warning_message[ERROR_LENGTH];
    sprintf (warning_message,"Attempting to open an open file %s",
	     name_.c_str());
    WARNING("FileHdf5::open",warning_message);

  } else {

    std::string full_name = path_ + "/" + name_;

    if (mode_ == "r") {

      file_id_ = H5Fopen(full_name.c_str(), 
		      H5F_ACC_RDONLY, 
		      H5P_DEFAULT);

    } else if (mode_ == "w") {

      file_id_ = H5Fcreate(full_name.c_str(), 
			H5F_ACC_TRUNC, 
			H5P_DEFAULT, 
			H5P_DEFAULT);

    } else {

      char error_message[ERROR_LENGTH];
      sprintf (error_message,"Unrecognized mode: %s",
	       mode_.c_str());
      ERROR("FileHdf5::open",error_message);

    }

    if (file_id_ >= 0) {

      is_file_open_ = true;

    } else {

      char warning_message[ERROR_LENGTH];
      sprintf (warning_message,
	       "Return value %d opening file %s",
	       file_id_,full_name.c_str());
      WARNING("FileHdf5::open",warning_message);

    }
  }
}

//----------------------------------------------------------------------

void FileHdf5::close () throw()
/**
 */
{
  std::string full_name = path_ + "/" + name_;

  if (! is_file_open_) {
    char warning_message[ERROR_LENGTH];
    sprintf (warning_message,
	     "Attempting to close a closed file %s",
	     full_name.c_str());
    WARNING("FileHdf5::close",warning_message);
  } else {
    int retval = H5Fclose (file_id_);
    if (retval >= 0) {
      is_file_open_ = false;
    } else {
      char warning_message[ERROR_LENGTH];
      sprintf (warning_message,
	       "Return value %d closing file %s",
	       retval,full_name.c_str());
      WARNING("FileHdf5::close",warning_message);
    }
  }
}

// //----------------------------------------------------------------------

// void FileHdf5::open_group (std::string group) throw()
// /**
//  */
// {
//   INCOMPLETE("FileHdf5::open_group");
// }

// //----------------------------------------------------------------------

// void FileHdf5::close_group () throw()
// /**
//  */
// {
//   INCOMPLETE("FileHdf5::open_group");
// }

//----------------------------------------------------------------------

void FileHdf5::data_set
(
 std::string      name,
 enum scalar_type type,
 int n0,  int n1,  int n2, int n3, int n4
 ) throw()
{
  if (is_data_open_) data_close_();

  if (mode_ != "w") {

    char error_message[ERROR_LENGTH];
    sprintf (error_message, "Expecting data_set() call for reading");
    ERROR("FileHdf5::data_set",error_message);

  } else {

    // Determine the dataset rank 
    // ASSUMES ( n[j] == 0 ) implies ( n[i] == 0 for i<j )

    data_name_ = name;

    data_rank_ = 5;

    if (n4 == 0) -- data_rank_;
    if (n3 == 0) -- data_rank_;
    if (n2 == 0) -- data_rank_;
    if (n1 == 0) -- data_rank_;

    ASSERT ("FileHdf5::data_set","d is out of range",
	    (1 <= data_rank_ && data_rank_ <= 5));

    if (n0) data_size_[0] = n0;
    if (n1) data_size_[1] = n1;
    if (n2) data_size_[2] = n2;
    if (n3) data_size_[3] = n3;
    if (n4) data_size_[4] = n4;

    // Define the dataspace
    data_space_id_ = H5Screate_simple (data_rank_, data_size_, NULL);
    
    // Create the dataset
    data_set_id_ = H5Dcreate( file_id_, 
			      data_name_.c_str(), 
			      type_(type), 
			      data_space_id_,  
			      H5P_DEFAULT );

    if (data_set_id_ < 0) {

     
    } else {

      char warning_message[ERROR_LENGTH];
      sprintf (warning_message,
	       "Return value %d opening dataset %s",
	       data_set_id_,data_name_.c_str());
      WARNING("FileHdf5::data_set",warning_message);
    }
  }

  data_open_();
}

//----------------------------------------------------------------------

void FileHdf5::data_get
(
 std::string name, 
 enum scalar_type * type, 
 int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  // is file open?
  // is data open?
  // return data size
}

//----------------------------------------------------------------------

void FileHdf5::data_read (void * buffer) throw()
{
  if (! is_file_open_) {
      char error_message[ERROR_LENGTH];
      sprintf (error_message,
	       "Trying to read from unopened file %s",
	       (path_ + name_).c_str());
      ERROR("FileHdf5::data_read",error_message);
  }
  if (! is_data_open_) {
  }
  // is file open?
  // is data open?
  // read data
}

//----------------------------------------------------------------------

void FileHdf5::data_write (const void * buffer) throw()
{
  // is file open?
  // is data open?
  // write data
}

//======================================================================

int FileHdf5::type_(enum scalar_type type) throw()
{
  // (*) NATIVE    -   FLOAT DOUBLE LDOUBLE
  // ( ) IEEE      -   F32BE F64BE     -
  // ( ) STD     B16BE B32BE B64BE     -
  // Types: http://www.hdfgroup.org/HDF5/Tutor/datatypes.html#native-types
  // char          H5T_NATIVE_CHAR   H5T_STD_I8BE or H5T_STD_I8LE
  // float         H5T_NATIVE_FLOAT   H5T_IEEE_F32BE or H5T_IEEE_F32LE  
  // double        H5T_NATIVE_DOUBLE   H5T_IEEE_F64BE or H5T_IEEE_F64LE  
  // unsigned char H5T_NATIVE_UCHAR   H5T_STD_U8BE or H5T_STD_U8LE
  // int           H5T_NATIVE_INT   H5T_STD_I32BE or H5T_STD_I32LE
  // short:        H5T_NATIVE_SHORT   H5T_STD_I16BE or H5T_STD_I16LE
  // long:         H5T_NATIVE_LONG   H5T_STD_I32BE, H5T_STD_I32LE,
  //               H5T_STD_I64BE or H5T_STD_I64LE
  // long long:    H5T_NATIVE_LLONG   H5T_STD_I64BE or H5T_STD_I64LE

  switch (type) {
  case scalar_type_unknown:
    ERROR("FileHdf5::type_",
	  "scalar_type_unknown not implemented");
    return 0;
    break;
  case scalar_type_float:
    return H5T_NATIVE_FLOAT;
    break;
  case scalar_type_double:
    return H5T_NATIVE_DOUBLE;
    break;
  case scalar_type_long_double:
    return H5T_NATIVE_LDOUBLE;
    break;
  case scalar_type_char:
    return H5T_NATIVE_CHAR;
    break;
  scalar_type_int:
    return H5T_NATIVE_INT;
    break;
  scalar_type_long_int:
    return H5T_NATIVE_LINT;
    break;
  }
}

//----------------------------------------------------------------------

void FileHdf5::data_open_() throw()
{
}

//----------------------------------------------------------------------

void FileHdf5::data_close_() throw()
{
}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void FileHdf5::open_data
(
 std::string data, 
 int *       nx,
 int *       ny,
 int *       nz
 ) throw()
/**
 */
{
  if (mode_ != "r") {

    char error_message[ERROR_LENGTH];
    sprintf (error_message, "Expecting data_set_id_open_write()");
    ERROR("FileHdf5::open_dataset",error_message);

  } else {

    int d;        // Dataset dimension
    hsize_t n[3]; // Dataset size

    if (mode_ == "r") {

      // Open the dataset for reading

      // Open the dataset
      data_set_id_ = H5Dopen( file_id_, data.c_str());
      // Get the size of dataset
      data_space_id_ = H5Dget_space (data_set_id_);
      d = H5Sget_simple_extent_ndims(data_space_id_);
      if (d > 3) {
	char error_message[ERROR_LENGTH];
	sprintf (error_message, "Dataset has too many dimensions %d",d);
	ERROR("FileHdf5::open_dataset",error_message);
      }

      H5Sget_simple_extent_dims(data_space_id_,n,0);

      // Set the array size accordingly
    
      *nx = n[0];
      *ny = (d > 1) ? n[1] : 1;
      *nz = (d > 2) ? n[2] : 1;

    }

    if (data_set_id_ >= 0) {

      data_name_ = data;
      
    } else {

      char warning_message[ERROR_LENGTH];
      sprintf (warning_message,
	       "Return value %d opening dataset %s",
	       data_set_id_,data.c_str());
      WARNING("FileHdf5::open_dataset",warning_message);
    }
  }
}

//----------------------------------------------------------------------

// void FileHdf5::dataset_set_hyperslab ()
// {
//   hsize_t offset[3] = {ix0,iy0,iz0};
//   hsize_t count[3]  = {mx,my,mz};
//   status_id_ = H5Sselect_hyperslab 
//     (data_space_id_, H5S_SELECT_SET,  offset, NULL, count, NULL);
// }

//----------------------------------------------------------------------

void FileHdf5::close_data () throw()
{ H5Dclose( data_set_id_); }

//----------------------------------------------------------------------

void FileHdf5::read  (char              * buffer,
		      enum precision_enum precision) throw()
/// @param buffer     Pointer to data buffer
/// @param precision  Precision of values in data buffer
{
  H5Dread (data_set_id_, 
	   type_(precision), 
	   data_space_id_, 
	   H5S_ALL, 
	   H5P_DEFAULT, 
	   buffer);
}

//----------------------------------------------------------------------

void FileHdf5::write (const char        * buffer,
		      enum precision_enum precision) throw()
/**
 */
{
  H5Dwrite (data_set_id_, 
	    type_(precision), 
	    data_space_id_, 
	    H5S_ALL, 
	    H5P_DEFAULT, 
	    buffer);
}


