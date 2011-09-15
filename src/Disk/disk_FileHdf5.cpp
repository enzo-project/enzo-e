// See LICENSE_CELLO file for license and copyright information

/// @file      disk_FileHdf5.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:11:36 PST 2008
/// @todo      Factor out common code from file_read_meta() and data_read_meta()
/// @todo      Factor out common code from file_write_meta() and data_write_meta()
/// @todo      Factor out common code from file_read_meta() and data_open(), data_read()
/// @todo      Factor out common code from file_write_meta() and data_create(), data_write()
/// @brief     Implementation of the FileHdf5 class

#include "cello.hpp"

#include "disk.hpp"

//----------------------------------------------------------------------
 
FileHdf5::FileHdf5 (std::string path, std::string name) throw()
  : File(path,name),
    file_id_(0),
    is_file_open_(false),
    data_set_id_(0),
    data_space_id_(0),
    attribute_id_(0),
    group_id_(0),
    group_name_("/"),
    is_group_open_(false),
    data_name_(""),
    data_type_(scalar_type_unknown),
    data_rank_(0),
    is_data_open_(false)
{
  const int rank_max = 5;
  for (int i=0; i<rank_max; i++) {
    data_size_[i] = 0;
  }
}

//----------------------------------------------------------------------

void FileHdf5::file_open () throw()
{
  if (is_file_open_) {

    WARNING1("FileHdf5::file_open",
	     "Attempting to reopen an opened file %s",
	     name_.c_str());

  } else {

    std::string full_name = path_ + "/" + name_;

    file_id_ = H5Fopen(full_name.c_str(), 
		       H5F_ACC_RDONLY, 
		       H5P_DEFAULT);

    // error check

    if (file_id_ >= 0) {

      is_file_open_ = true;

    } else {

      WARNING2("FileHdf5::file_open",
	       "Return value %d opening file %s",
	       file_id_,full_name.c_str());

    }
  }
}

//----------------------------------------------------------------------

void FileHdf5::file_create () throw()
{

  std::string full_name = path_ + "/" + name_;

  file_id_ = H5Fcreate(full_name.c_str(), 
		       H5F_ACC_TRUNC, 
		       H5P_DEFAULT, 
		       H5P_DEFAULT);

  // error check

  if (file_id_ >= 0) {

    is_file_open_ = true;

  } else {

    WARNING2("FileHdf5::file_create",
	    "Return value %d opening file %s",
	    file_id_,full_name.c_str());

  }
}

//----------------------------------------------------------------------

void FileHdf5::file_close () throw()
{

  int retval;

  if (is_data_open_) {

    // Close dataspace

    retval = H5Sclose (data_space_id_);

    // error check

    if (retval >= 0) {

      is_data_open_ = false;

    } else {

      WARNING2("FileHdf5::file_close",
	      "Return value %d closing dataspace %s",
	      retval,data_name_.c_str());

    }

    // Close dataset

    retval = H5Dclose (data_set_id_);

    // error check

    if (retval >= 0) {

      is_data_open_ = false;

    } else {

      WARNING2("FileHdf5::file_close",
	      "Return value %d closing dataset %s",
	      retval,data_name_.c_str());
    }

  }

  if (is_file_open_) {

    int retval = H5Fclose (file_id_);

    if (retval >= 0) {

      is_file_open_ = false;

    } else {

      WARNING2("FileHdf5::file_close",
	      "Return value %d closing file %s",
	       retval,(path_ + "/" + name_).c_str());
    }

  } else {

    // error check

    WARNING1("FileHdf5::file_close",
	    "Closing an already closed file %s",
	     (path_ + "/" + name_).c_str());
  }
}

//----------------------------------------------------------------------

void FileHdf5::data_open
( std::string name,  enum scalar_type * type,
  int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{

 // Error checking

  if (! is_file_open_) {

    ERROR1("FileHdf5::data_read",
	  "Trying to read from unopened file %s",
	  (path_ + name_).c_str());
  }

  // Open the data set
  
  hid_t group = (is_group_open_) ? group_id_ : file_id_;

  data_set_id_ = H5Dopen( group, name.c_str());

  // error check

  if (data_set_id_ < 0) {

    ERROR3("FileHdf5::data_open",
	   "H5Dopen() returned %d when opening dataset %s in file %s",
	   data_set_id_,name.c_str(),(path_ + name_).c_str());
  }

  is_data_open_ = true;

  // Get dataset size

  data_space_id_ = H5Dget_space (data_set_id_);

  // error check rank

  int rank = H5Sget_simple_extent_ndims(data_space_id_);

  if (rank > 5) {

    ERROR3("FileHdf5::data_get",
	  "Dataset %s in file %s has unsupported rank %d",
	  name.c_str(),name_.c_str(),rank);
  }

  // get extents

  hsize_t n[5];
  H5Sget_simple_extent_dims(data_space_id_,n,0);

  // Initialize name and type

  data_name_ = name;
  data_type_ = hdf5_to_scalar_(H5Dget_type (data_set_id_));

  // set output parameters

  if (type) (*type) = data_type_;

  if (n0) (*n0) = n[0];
  if (n1) (*n1) = n[1];
  if (n2) (*n2) = n[2];
  if (n3) (*n3) = n[3];
  if (n4) (*n4) = n[4];

}

//----------------------------------------------------------------------

void FileHdf5::data_create
( std::string name,  enum scalar_type type,
  int n0, int n1, int n2, int n3, int n4) throw()
{

  // Initialize data attributes

  data_name_ = name;
  data_type_ = type;

  data_rank_ = 5;

  if (n4 == 0) -- data_rank_;
  if (n3 == 0) -- data_rank_;
  if (n2 == 0) -- data_rank_;
  if (n1 == 0) -- data_rank_;

  ASSERT ("FileHdf5::data_set","d is out of range",
	  (1 <= data_rank_ && data_rank_ <= 5));

  data_size_[0] = n0;
  data_size_[1] = n1;
  data_size_[2] = n2;
  data_size_[3] = n3;
  data_size_[4] = n4;

  // Define the dataspace                                                      

  data_space_id_ = H5Screate_simple (data_rank_, data_size_, NULL);

  // Create the new dataset

  hid_t group = (is_group_open_) ? group_id_ : file_id_;

  data_set_id_ = H5Dcreate( group,
			    name.c_str(),
			    scalar_to_hdf5_(type),
			    data_space_id_,
			    H5P_DEFAULT );

  // error check

  if (data_set_id_ < 0) {

    WARNING2("FileHdf5::data_set",
	    "Return value %d opening dataset %s",
	    data_set_id_,name.c_str());
  }

  // Update data state

  data_type_ = type;
  is_data_open_ = true;


}

//----------------------------------------------------------------------

void FileHdf5::data_read
( void * buffer) throw()
{

  // error check file open

  if (! is_file_open_) {

    ERROR1("FileHdf5::data_read",
	  "Trying to read from unopened file %s",
	  (path_ + name_).c_str());
  }

  // error check dataset open

  if (! is_data_open_) {

    ERROR1("FileHdf5::data_read",
	  "Trying to read unopened dataset %s",
	  data_name_.c_str());
  }

  // read data

  H5Dread (data_set_id_,
           scalar_to_hdf5_(data_type_),
           data_space_id_,
           H5S_ALL,
           H5P_DEFAULT,
           buffer);

}

//----------------------------------------------------------------------

void FileHdf5::data_write
( const void * buffer) throw()
{
  // error check file open

  if (! is_file_open_) {

    ERROR1("FileHdf5::data_write",
	  "Trying to write to unopened file %s",
	  (path_ + name_).c_str());
  }

  // error check dataset open

  if (! is_data_open_) {

    ERROR1("FileHdf5::data_write",
	  "Trying to write unopened dataset %s",
	  data_name_.c_str());
  }

  H5Dwrite (data_set_id_,
            scalar_to_hdf5_(data_type_),
            data_space_id_,
            H5S_ALL,
            H5P_DEFAULT,
            buffer);
}

//----------------------------------------------------------------------

void FileHdf5::data_close() throw()
{
  if (is_data_open_) {
    H5Sclose (data_space_id_);
    H5Dclose( data_set_id_);
    is_data_open_ = false;
  }
}

//----------------------------------------------------------------------

void FileHdf5::file_read_meta
  ( void * buffer, std::string name,  enum scalar_type * type,
    int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  if ( ! is_file_open_) {

    ERROR1("FileHdf5::file_read_meta",
	  "Trying to read metadata from the unopened file %s",
	  (path_ + name_).c_str());
  }

  hid_t meta_id = H5Aopen_name(file_id_, name.c_str());

  if (meta_id < 0) {

    ERROR3("FileHdf5::file_read_meta",
	  "H5Aopen_name() returned %d when opening attribute %s in file %s",
	  meta_id, name.c_str(),(path_ + name_).c_str());
  }

  // Get attribute size

  hid_t meta_space_id = H5Aget_space (meta_id);

  // error check rank

  int rank = H5Sget_simple_extent_ndims(meta_space_id);

  if (rank > 5) {
    ERROR3("FileHdf5::file_read_meta",
	  "Attribute %s in file %s has unsupported rank %d",
	  name.c_str(),name_.c_str(),rank);
  }

  // set output parameters

  scalar_type scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

  if (type) (*type) = scalar_type;

  hsize_t n[5];
  H5Sget_simple_extent_dims(meta_space_id,n,0);

  if (n0) (*n0) = n[0];
  if (n1) (*n1) = n[1];
  if (n2) (*n2) = n[2];
  if (n3) (*n3) = n[3];
  if (n4) (*n4) = n[4];

  // Read the attribute

  H5Aread(meta_id, scalar_to_hdf5_(scalar_type), buffer);

}

//----------------------------------------------------------------------

void FileHdf5::file_write_meta
  ( const void * buffer, std::string name, enum scalar_type type,
    int n0, int n1, int n2, int n3, int n4) throw()
{
  if ( ! is_file_open_) {
    ERROR1("FileHdf5::file_write_meta",
	  "Trying to write metadata to the unopened file %s",
	  (path_ + name_).c_str());
  }

  // Determine the attribute rank

  hsize_t meta_rank = 5;

  if (n4 == 0) -- meta_rank;
  if (n3 == 0) -- meta_rank;
  if (n2 == 0) -- meta_rank;
  if (n1 == 0) -- meta_rank;

  ASSERT ("FileHdf5::data_set","d is out of range",
	  (1 <= meta_rank && meta_rank <= 5));


  hsize_t meta_size[5];

  meta_size[0] = n0;
  meta_size[1] = n1;
  meta_size[2] = n2;
  meta_size[3] = n3;
  meta_size[4] = n4;
  
  // Open the data space

  hid_t meta_space_id = H5Screate_simple (meta_rank, meta_size, NULL);

  hid_t meta_id = H5Acreate ( file_id_,
			      name.c_str(),
			      scalar_to_hdf5_(type),
			      meta_space_id,
			      H5P_DEFAULT);

  H5Awrite (meta_id, scalar_to_hdf5_(type), buffer);

  H5Sclose(meta_space_id);
  H5Aclose(meta_id);
}

//----------------------------------------------------------------------

void FileHdf5::data_read_meta
  ( void * buffer, std::string name,  enum scalar_type * type,
    int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  // error check file open

  if ( ! is_file_open_) {

    ERROR1("FileHdf5::data_read_meta",
	  "Trying to read attribute from the unopened file %s",
	  (path_ + name_).c_str());
  }

  // error check dataset open

  if (! is_data_open_) {

    ERROR1("FileHdf5::data_read_meta",
	  "Trying to read attribute from unopened dataset %s",
	  data_name_.c_str());
  }

  hid_t meta_id = H5Aopen_name(data_set_id_, name.c_str());

  if (meta_id < 0) {

    ERROR3("FileHdf5::data_read_meta",
	  "H5Aopen_name() returned %d when opening attribute %s in file %s",
	  meta_id, name.c_str(),(path_ + name_).c_str());
  }

  // Get attribute size

  hid_t meta_space_id = H5Aget_space (meta_id);

  // error check rank

  int rank = H5Sget_simple_extent_ndims(meta_space_id);

  if (rank > 5) {

    ERROR3("FileHdf5::data_read_meta",
	  "Attribute %s in file %s has unsupported rank %d",
	  name.c_str(),name_.c_str(),rank);
  }

  // set output parameters

  scalar_type scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

  if (type) (*type) = scalar_type;

  hsize_t n[5];
  H5Sget_simple_extent_dims(meta_space_id,n,0);

  if (n0) (*n0) = n[0];
  if (n1) (*n1) = n[1];
  if (n2) (*n2) = n[2];
  if (n3) (*n3) = n[3];
  if (n4) (*n4) = n[4];

  // Read the attribute

  H5Aread(meta_id, scalar_to_hdf5_(scalar_type), buffer);

}

//----------------------------------------------------------------------

void FileHdf5::data_write_meta
  ( const void * buffer, std::string name, enum scalar_type type,
    int n0, int n1, int n2, int n3, int n4) throw()
{
  // error check file open

  if ( ! is_file_open_) {

    ERROR1("FileHdf5::data_write_meta",
	  "Trying to write metadata to the unopened file %s",
	  (path_ + name_).c_str());
  }

  // error check dataset open

  if (! is_data_open_) {

    ERROR1("FileHdf5::data_write_meta",
	  "Trying to read attribute from unopened dataset %s",
	  data_name_.c_str());
  }

  // Determine the attribute rank

  hsize_t meta_rank = 5;

  if (n4 == 0) -- meta_rank;
  if (n3 == 0) -- meta_rank;
  if (n2 == 0) -- meta_rank;
  if (n1 == 0) -- meta_rank;

  ASSERT ("FileHdf5::data_set","d is out of range",
	  (1 <= meta_rank && meta_rank <= 5));


  hsize_t meta_size[5];

  meta_size[0] = n0;
  meta_size[1] = n1;
  meta_size[2] = n2;
  meta_size[3] = n3;
  meta_size[4] = n4;
  
  // Open the data space

  hid_t meta_space_id = H5Screate_simple (meta_rank, meta_size, NULL);

  hid_t meta_id = H5Acreate ( data_set_id_,
			      name.c_str(),
			      scalar_to_hdf5_(type),
			      meta_space_id,
			      H5P_DEFAULT);

  H5Awrite (meta_id, scalar_to_hdf5_(type), buffer);

  H5Sclose(meta_space_id);
  H5Aclose(meta_id);
}  

//----------------------------------------------------------------------

void FileHdf5::group_open (std::string name) throw()
{
  // Close current group if open
  group_close();
  
  // Error if not an absolute path name
  if (name[0] != '/') {

    ERROR1("FileHdf5::data_write_meta",
	  "Group name '%s' must begin with '/'", name.c_str());
  }
  
  group_id_ = H5Gopen(file_id_, name.c_str());

  is_group_open_ = true;
  
}

//----------------------------------------------------------------------

void FileHdf5::group_create (std::string group_path) throw()
{
  // Close current group if open
  group_close();

  // Error if not an absolute path name
  if (group_path[0] != '/') {

    ERROR1("FileHdf5::data_write_meta",
	  "Group name '%s' must begin with '/'", group_path.c_str());
  }


  // Start at root group
  
  std::string group_full = "/";
  std::string group_rest = group_path;
  group_rest.erase(0,1);

  group_id_ = H5Gopen(file_id_,group_full.c_str());

  // Loop through ancestor groups, creating if needed

  size_t pos;

  bool done = false;
  while ( ! done ) {

    pos = group_rest.find("/",0);

    // Get the next subgroup name

    std::string group = group_rest.substr(0,pos);
    group_rest.erase(0,pos+1);

    // Loop through children to find if subgroup exists

    H5G_info_t group_info;
    H5Gget_info_by_name (file_id_, group_full.c_str(), &group_info, H5P_DEFAULT);

    bool group_exists = false;
    for (int i=0; i<group_info.nlinks; i++) {
      char group_name[80];
      H5Lget_name_by_idx (file_id_, group_full.c_str(),
			  H5_INDEX_NAME,
			  H5_ITER_NATIVE,i,
			  group_name,80,H5P_DEFAULT);
      if (group == group_name) {
	group_exists = true;
	break;
      }
    }

    group_full = group_full + group + "/" ;

    //  Open or create next group in path

    hid_t group_new;

    if (group_exists) {
      group_new= H5Gopen(file_id_,group_full.c_str());
    } else {
      group_new = H5Gcreate (file_id_,group_full.c_str(), H5P_DEFAULT);
    }

    // Close parent group

    H5Gclose (group_id_);

    // Update group id

    group_id_ = group_new;
    
    done = (pos == std::string::npos);

  }

  is_group_open_ = true;
}

//----------------------------------------------------------------------

void FileHdf5::group_close () throw()
{
  herr_t status;
  if (is_group_open_) {
    //
    status = H5Gclose(group_id_);
  }
  is_group_open_ = false;
}

//======================================================================

int FileHdf5::scalar_to_hdf5_ (enum scalar_type type) const throw()
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
    ERROR1("FileHdf5::type_", "unsupported type %d", type);
    hdf5_type = 0;
    break;
  }
  return hdf5_type;
}

//----------------------------------------------------------------------

enum scalar_type FileHdf5::hdf5_to_scalar_ (int hdf5_type) const throw()
{

  H5T_class_t hdf5_class = H5Tget_class(hdf5_type);
  size_t      hdf5_size  = H5Tget_size (hdf5_type);

  enum scalar_type type;

  if (hdf5_class == H5T_INTEGER) {

    if (hdf5_size == sizeof(char))   type = scalar_type_char;
    if (hdf5_size == sizeof(int))    type = scalar_type_int;
    if (hdf5_size == sizeof(long))   type = scalar_type_long;

  } else if (hdf5_class == H5T_FLOAT) {

    if (hdf5_size == sizeof(float))  type = scalar_type_float;
    if (hdf5_size == sizeof(double)) type = scalar_type_double;

  } else {

    ERROR2("FileHdf5::data_get",
	  "Unknown type of class %d and size %d",
	  hdf5_class, int(hdf5_size));
  }

  return type;

}

