// See LICENSE_CELLO file for license and copyright information

/// @file      disk_FileHdf5.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:11:36 PST 2008
/// @brief     Implementation of the FileHdf5 class

#include "cello.hpp"

#include "disk.hpp"
 
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
  int retval;
  if (is_data_open_) {
    retval = H5Sclose (data_space_id_);
    if (retval >= 0) {
      is_data_open_ = false;
    } else {
      char warning_message[ERROR_LENGTH];
      sprintf (warning_message,
	       "Return value %d closing dataspace %s",
	       retval,data_name_.c_str());
      WARNING("FileHdf5::close",warning_message);
    }
    retval = H5Dclose (data_set_id_);
    if (retval >= 0) {
      is_data_open_ = false;
    } else {
      char warning_message[ERROR_LENGTH];
      sprintf (warning_message,
	       "Return value %d closing dataset %s",
	       retval,data_name_.c_str());
      WARNING("FileHdf5::close",warning_message);
    }
  }
  if (is_file_open_) {
    int retval = H5Fclose (file_id_);
    if (retval >= 0) {
      is_file_open_ = false;
    } else {
      std::string full_name = path_ + "/" + name_;
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

void FileHdf5::data_get_
(
 std::string name, 
 enum scalar_type * type, 
 int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  // Error checking

  if (! is_file_open_) {
    char error_message[ERROR_LENGTH];
    sprintf (error_message,
	     "Trying to read from unopened file %s",
	     (path_ + name_).c_str());
    ERROR("FileHdf5::data_read",error_message);
  }

  data_name_ = name;

  // Open the data set
  
  data_set_id_ = H5Dopen( file_id_, name.c_str());
  is_data_open_ = true;

  // Get dataset size

  data_space_id_ = H5Dget_space (data_set_id_);

  // get rank

  int rank = H5Sget_simple_extent_ndims(data_space_id_);

  if (rank > 5) {
    char error_message[ERROR_LENGTH];
    sprintf (error_message, 
	     "Dataset %s in file %s has unsupported rank %d",
	     name.c_str(),name_.c_str(),rank);
    ERROR("FileHdf5::data_get",error_message);
  }

  // get extents

  hsize_t n[5];
  H5Sget_simple_extent_dims(data_space_id_,n,0);

  // set output parameters

  if (type) {

    hid_t hdf5_type = H5Dget_type (data_set_id_);

    H5T_class_t hdf5_class = H5Tget_class(hdf5_type);
    size_t      hdf5_size  = H5Tget_size(hdf5_type);

    if (hdf5_class == H5T_INTEGER && hdf5_size == sizeof(char)) {
      (*type) = scalar_type_char;
    } else if (hdf5_class == H5T_INTEGER && hdf5_size == sizeof(int)) {
      (*type) = scalar_type_int;
    } else if (hdf5_class == H5T_INTEGER && hdf5_size == sizeof(long)) {
      (*type) = scalar_type_long;
    } else if (hdf5_class == H5T_FLOAT && hdf5_size == sizeof(float)) {
      (*type) = scalar_type_float;
    } else if (hdf5_class == H5T_FLOAT && hdf5_size == sizeof(double)) {
      (*type) = scalar_type_double;
    } else {
      char error_message[ERROR_LENGTH];
      sprintf (error_message, 
	       "Unknown type of class %d and size %d",
	       hdf5_class, hdf5_size);
      ERROR("FileHdf5::data_get",error_message);
    }
  }

  data_type_ = (*type);

  if (n0) (*n0) = n[0];
  if (n1) (*n1) = n[1];
  if (n2) (*n2) = n[2];
  if (n3) (*n3) = n[3];
  if (n4) (*n4) = n[4];

}

//----------------------------------------------------------------------

void FileHdf5::data_set_
(
 std::string      name,
 enum scalar_type type,
 int n0,  int n1,  int n2, int n3, int n4
 ) throw()
{
  // Close the previous data set if its open
  if (is_data_open_) {
    H5Sclose (data_space_id_);
    H5Dclose( data_set_id_);
    is_data_open_ = false;
  }

  if (mode_ != "w") {

    char error_message[ERROR_LENGTH];
    sprintf (error_message, "Expecting data_set() call for reading");
    ERROR("FileHdf5::data_set",error_message);

  } else {

    // Determine the dataset rank 
    // ASSUMES ( n[j] == 0 ) implies ( n[i] == 0 for i<j )

    data_name_ = name;
    data_type_ = type;

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
			      name.c_str(), 
			      hdf5_type_(type), 
			      data_space_id_,  
			      H5P_DEFAULT );

    if (data_set_id_ < 0) {

      char warning_message[ERROR_LENGTH];
      sprintf (warning_message,
	       "Return value %d opening dataset %s",
	       data_set_id_,name.c_str());
      WARNING("FileHdf5::data_set",warning_message);
    }
  }

  // Open the data set

  // data_set_id_ = H5Dopen( file_id_, name.c_str());

  data_type_ = type;
  is_data_open_ = true;
}

//----------------------------------------------------------------------
void FileHdf5::data_read
(
 void * buffer,
 std::string name, 
 enum scalar_type * type, 
 int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  data_get_  (name,type,n0,n1,n2,n3,n4);
  data_read_ (buffer);
}

//----------------------------------------------------------------------

void FileHdf5::data_read_ (void * buffer) throw()
{
  if (! is_file_open_) {
    char error_message[ERROR_LENGTH];
    sprintf (error_message,
	     "Trying to read from unopened file %s",
	     (path_ + name_).c_str());
    ERROR("FileHdf5::data_read",error_message);
  }

  if (! is_data_open_) {
    char error_message[ERROR_LENGTH];
    sprintf (error_message,
	     "Trying to read unopened dataset %s",
	     data_name_.c_str());
    ERROR("FileHdf5::data_read",error_message);
  }

  H5Dread (data_set_id_, 
 	   hdf5_type_(data_type_), 
 	   data_space_id_, 
 	   H5S_ALL, 
 	   H5P_DEFAULT, 
 	   buffer);
}

//----------------------------------------------------------------------

void FileHdf5::data_write
(
 const void * buffer,
 std::string      name,
 enum scalar_type type,
 int n0,  int n1,  int n2, int n3, int n4
 ) throw()
{
  data_set_   (name,type,n0,n1,n2,n3,n4);
  data_write_ (buffer);
}

//----------------------------------------------------------------------

void FileHdf5::data_write_ (const void * buffer) throw()
{
  if (! is_file_open_) {
    char error_message[ERROR_LENGTH];
    sprintf (error_message,
	     "Trying to write to unopened file %s",
	     (path_ + name_).c_str());
    ERROR("FileHdf5::data_write",error_message);
  }

  if (! is_data_open_) {
    char error_message[ERROR_LENGTH];
    sprintf (error_message,
	     "Trying to write unopened dataset %s",
	     data_name_.c_str());
    ERROR("FileHdf5::data_write",error_message);
  }

  // is file open?
  // is data open?
  // write data
  H5Dwrite (data_set_id_, 
   	    hdf5_type_(data_type_), 
   	    data_space_id_, 
   	    H5S_ALL, 
   	    H5P_DEFAULT, 
   	    buffer);

  // hid_t dataset_id, 
  //   hid_t mem_type_id, 
  //   hid_t mem_space_id, 
  //   hid_t file_space_id, 
  //   hid_t xfer_plist_id, 
  //   const void * buf ) 
}

//======================================================================

int FileHdf5::hdf5_type_(enum scalar_type type) throw()
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

  hid_t hdf5_type;
  switch (type) {
  case scalar_type_unknown:
    ERROR("FileHdf5::type_",
	  "scalar_type_unknown not implemented");
    hdf5_type = 0;
    break;
  case scalar_type_float:
    hdf5_type = H5T_NATIVE_FLOAT;
    //  H5T_IEEE_F32BE 
    //  H5T_IEEE_F32LE  
    break;
  case scalar_type_double:
    hdf5_type = H5T_NATIVE_DOUBLE;
    //  H5T_IEEE_F64BE
    //  H5T_IEEE_F64LE  
    break;
  case scalar_type_long_double:
    ERROR("FileHdf5::type_",
	  "long double not supported");
    hdf5_type = 0;
    break;
  case scalar_type_char:
    hdf5_type = H5T_NATIVE_CHAR;
    //  H5T_STD_I8BE
    //  H5T_STD_I8LE
    break;
  case scalar_type_int:
    hdf5_type = H5T_NATIVE_INT;
    //  H5T_STD_I32BE 
    //  H5T_STD_I32LE
    break;
  case scalar_type_long:
    //  H5T_STD_I64BE
    //  H5T_STD_I64LE
    hdf5_type = H5T_NATIVE_LONG;
    break;
  default:
    break;
  }
  return hdf5_type;
}

//----------------------------------------------------------------------

// void FileHdf5::dataset_set_hyperslab ()
// {
//   hsize_t offset[3] = {ix0,iy0,iz0};
//   hsize_t count[3]  = {mx,my,mz};
//   status_id_ = H5Sselect_hyperslab 
//     (data_space_id_, H5S_SELECT_SET,  offset, NULL, count, NULL);
// }



