// See LICENSE_CELLO file for license and copyright information

/// @file      disk_FileHdf5.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:11:36 PST 2008
/// @brief     Implementation of the FileHdf5 class

#include "cello.hpp"

#include "disk.hpp"

#define MAX_DATA_RANK 3
#define MAX_ATTR_RANK 3

//----------------------------------------------------------------------
 
FileHdf5::FileHdf5 (std::string path, std::string name) throw()
  : File(path,name),
    file_id_(0),
    is_file_open_(false),
    data_id_(0),
    data_space_id_(H5S_ALL),
    mem_space_id_(H5S_ALL),
    attribute_id_(0),
    group_id_(0),
    group_name_("/"),
    group_prop_(H5P_DEFAULT),
    is_group_open_(false),
    data_name_(""),
    data_type_(scalar_type_unknown),
    data_rank_(0),
    data_prop_(H5P_DEFAULT),
    is_data_open_(false),
    compress_level_(0)
{
  for (int i=0; i<MAX_DATA_RANK; i++) {
    data_dims_[i] = 0;
  }

  // group_prop_ = H5Pcreate (H5P_GROUP_CREATE);
  data_prop_  = H5Pcreate (H5P_DATASET_CREATE);
  group_prop_ = H5P_DEFAULT;
  //data_prop_ = H5P_DEFAULT;
}

//----------------------------------------------------------------------

FileHdf5::~FileHdf5() throw()
{
  // H5Pclose (group_prop_);
  H5Pclose (data_prop_);
}

//----------------------------------------------------------------------

void FileHdf5::file_open () throw()
{

  // check file closed

  std::string file_name = path_ + "/" + name_;

  ASSERT1("FileHdf5::file_open", "Attempting to reopen an opened file %s",
	  file_name.c_str(), ! is_file_open_);

  // open file

  file_id_ = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // error check file opened

  ASSERT2("FileHdf5::file_open", "Return value %d opening file %s",
	 file_id_,file_name.c_str(), file_id_ >= 0);

  // update file state

  is_file_open_ = true;

}

//----------------------------------------------------------------------

void FileHdf5::file_create () throw()
{

  // create file

  std::string file_name = path_ + "/" + name_;

  file_id_ = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // error check file created

  TRACE1("File = %s",file_name.c_str());
  ASSERT2("FileHdf5::file_create",  "Return value %d opening file %s",
	  file_id_,file_name.c_str(), file_id_ >= 0);

  // update file state

  is_file_open_ = true;

}

//----------------------------------------------------------------------

void FileHdf5::file_close () throw()
{
  // error check file open

  std::string file_name = path_ + "/" + name_;

  ASSERT1("FileHdf5::file_close", "Closing an already closed file %s",
	  file_name.c_str(), (is_file_open_));

  // close dataset if opened

  int retval;

  if (is_data_open_) {

    // close dataspace

    close_space_ (data_space_id_);
    close_space_ (mem_space_id_);

    // Close dataset

    close_dataset_ ();

    // update dataset state

    is_data_open_ = false;

  }
  // Close the file

  retval = H5Fclose (file_id_);

  // error check H5Fclose

  ASSERT2("FileHdf5::file_close", "Return value %d closing file %s",
	  retval,file_name.c_str(), (retval >= 0));

  // update file state

  is_file_open_ = false;
}

//----------------------------------------------------------------------

void FileHdf5::data_open
( std::string name,  scalar_type * type,
  int * nx, int * ny, int * nz) throw()
{

 // error check file closed

  std::string file_name = path_ + "/" + name_;

  ASSERT1("FileHdf5::data_open", "Trying to read from unopened file %s",
	  file_name.c_str(), is_file_open_ );

  // Open the dataset
  
  hid_t group = (is_group_open_) ? group_id_ : file_id_;

  data_id_ = open_dataset_(group,name);

  is_data_open_ = true;

  // Get the dataspace

  data_space_id_ = get_data_space_(data_id_, name);

  // set output extents

  get_output_extents_(data_space_id_,nx,ny,nz);

  // Initialize name and type

  data_name_ = name;
  data_type_ = hdf5_to_scalar_(H5Dget_type (data_id_));
  WARNING("FileHdf5::data_open","Set memspace");

  // set output parameters

  if (type) (*type) = data_type_;

}

//----------------------------------------------------------------------

void FileHdf5::data_create
( std::string name,  scalar_type type,
  int nxd, int nyd, int nzd,
  int nx, int ny, int nz) throw()
{

  if (nx==0) nx=nxd;
  if (ny==0) ny=nyd;
  if (nz==0) nz=nzd;

  // Initialize data attributes

  data_name_ = name;
  data_type_ = type;

  // Determine the group

  hid_t group = (is_group_open_) ? group_id_ : file_id_;

  // Create dataspace


  data_space_id_ = create_data_space_ (nxd,nyd,nzd,nx,ny,nz);
  mem_space_id_  = create_mem_space_  (nxd,nyd,nzd,nx,ny,nz);

  // Create the new dataset

  data_id_ = H5Dcreate( group,
			name.c_str(),
			scalar_to_hdf5_(type),
			data_space_id_,
			data_prop_ );

  // error check H5Dcreate

  ASSERT2("FileHdf5::data_create", "Return value %d creating dataset %s",
	  data_id_,name.c_str(), data_id_ >= 0);

  // update dataset state

  data_type_ = type;
  is_data_open_ = true;


}

//----------------------------------------------------------------------

void FileHdf5::data_read
( void * buffer) throw()
{

  // error check file open

  std::string file_name = path_ + "/" + name_;

  ASSERT1("FileHdf5::data_read", "Trying to read from unopened file %s",
	  file_name.c_str(), is_file_open_);

  // error check dataset open

  ASSERT1("FileHdf5::data_read", "Trying to read unopened dataset %s",
	  data_name_.c_str(), is_data_open_);

  // read data

  int retval = 
    H5Dread (data_id_,
	     scalar_to_hdf5_(data_type_),
	     data_space_id_,
	     mem_space_id_,
	     H5P_DEFAULT,
	     buffer);

  // error check H5Dread

  ASSERT1("FileHdf5::data_read","H5Dread() returned %d",retval,(retval>=0));

}

//----------------------------------------------------------------------

void FileHdf5::data_write ( const void * buffer ) throw()
{
  // error check file open

  std::string file_name = path_ + "/" + name_;

  ASSERT1("FileHdf5::data_write", "Trying to write to unopened file %s",
	  file_name.c_str(), is_file_open_);

  // error check dataset open

  ASSERT1("FileHdf5::data_write", "Trying to write unopened dataset %s",
	   data_name_.c_str(), (is_data_open_));

  // Write dataset to the file

  int retval = 
    H5Dwrite (data_id_,
	      scalar_to_hdf5_(data_type_),
	      mem_space_id_,
	      H5S_ALL,
	      H5P_DEFAULT,
	      buffer);

  // error check H5Dread

  ASSERT1("FileHdf5::data_write","H5Dwrite() returned %d",retval,(retval>=0));

}

//----------------------------------------------------------------------

void FileHdf5::data_close() throw()
{
  if (is_data_open_) {

    // close the dataspace

    close_space_(data_space_id_);
    close_space_(mem_space_id_);

    // close the dataset

    close_dataset_ ();

    // update dataset state

    is_data_open_ = false;
  }
}

//----------------------------------------------------------------------

void FileHdf5::file_read_meta
  ( void * buffer, std::string name,  scalar_type * type,
    int * nx, int * ny, int * nz) throw()
{

  std::string file_name = path_ + "/" + name_;

  // error check file open

  ASSERT1("FileHdf5::file_read_meta",
	  "Trying to read metadata from the unopened file %s",
	  file_name.c_str(), is_file_open_);

  hid_t meta_id = H5Aopen_name(file_id_, name.c_str());

  // error check H5Aopen_name

  ASSERT3("FileHdf5::file_read_meta",
	  "H5Aopen_name() returned %d opening %s in file %s",
	  meta_id, name.c_str(),file_name.c_str(),
	  (meta_id >= 0));

  // get dataspace

  hid_t meta_space_id = get_attr_space_(meta_id,name);

  // set output extents

  get_output_extents_(meta_space_id,nx,ny,nz);

  // set output type

  scalar_type scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

  if (type) (*type) = scalar_type;

  // Read the attribute

  int retval = 
    H5Aread(meta_id, scalar_to_hdf5_(scalar_type), buffer);

  // error check H5Aread

  ASSERT1("FileHdf5::file_read_meta_","H5Aread() returned %d",
	  retval,(retval>=0));

}

//----------------------------------------------------------------------

void FileHdf5::data_read_meta
  ( void * buffer, std::string name,  scalar_type * type,
    int * nx, int * ny, int * nz) throw()
{
  // error check file open

  std::string file_name = path_ + "/" + name_;

  ASSERT1("FileHdf5::data_read_meta",
	 "Trying to read attribute from the unopened file %s",
	  file_name.c_str(),
	  is_file_open_);

  // error check dataset open

  ASSERT1("FileHdf5::data_read_meta",
	  "Trying to read attribute from unopened dataset %s",
	  data_name_.c_str(),
	  is_data_open_);

  hid_t meta_id = H5Aopen_name(data_id_, name.c_str());

  ASSERT3("FileHdf5::data_read_meta",
	  "H5Aopen_name() returned %d when opening attribute %s in file %s",
	   meta_id, name.c_str(),file_name.c_str(),
	   (meta_id >= 0));

  // Get attribute size

  hid_t meta_space_id = get_attr_space_ (meta_id,name);

  // set output extents

  get_output_extents_(meta_space_id,nx,ny,nz);

  // set output parameters

  scalar_type scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

  if (type) (*type) = scalar_type;

  // Read the attribute

  int retval = H5Aread
    (meta_id, scalar_to_hdf5_(scalar_type), buffer);

  // error check H5Aread

  ASSERT1("FileHdf5::data_read_meta","H5Aread returned %d",retval,(retval>=0));
}

//----------------------------------------------------------------------

int FileHdf5::group_count () const throw()
{
  H5G_info_t group_info;
  H5Gget_info (group_id_, &group_info);
  return group_info.nlinks;
}

//----------------------------------------------------------------------

std::string FileHdf5::group_name (size_t i) const throw()
{
  char buffer[10];

  // 1.6.0 <= HDF5 version < 1.8.0
  //  H5Gget_objname_by_idx(group_id_,i,buffer,10);

  // 1.8.0 <= HDF5 version 
  H5Lget_name_by_idx (group_id_,group_name_.c_str(),H5_INDEX_NAME,H5_ITER_INC,
		      i,buffer,10,H5P_DEFAULT);

  return std::string(buffer);
}

//----------------------------------------------------------------------

void FileHdf5::group_chdir (std::string group_path) throw()
{
  // convert to absolute path if it is relative

  if (group_path[0] != '/') {
    group_path = relative_to_absolute_(group_path, group_name_);
  }

  // update the stored group name

  group_name_ = group_path;

  // error check group name

  ASSERT1("FileHdf5::group_chdir",
	  "Group name '%s' must begin with '/'", group_path.c_str(),
	  group_path[0] == '/');

}

//----------------------------------------------------------------------

void FileHdf5::group_open () throw()
{
  // close current group if open

  group_close();
  
  // open group

  group_id_ = H5Gopen(file_id_, group_name_.c_str());

  // error check H5Gopen()

  ASSERT2("FileHdf5::group_open()", "H5Gopen(%s) returned %d", 
	  group_id_,group_name_.c_str(),group_id_>=0);

  // update group state

  is_group_open_ = true;
  
}

//----------------------------------------------------------------------

void FileHdf5::group_create () throw()
{
  // close current group if open

  group_close();

  // Create ancestor groups beginning at root '/'
  
  std::string group_full = "/";
  std::string group_rest = group_name_;
  group_rest.erase(0,1);

  group_id_ = H5Gopen(file_id_,group_full.c_str());

  // loop through ancestor groups

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
    for (size_t i=0; i<group_info.nlinks; i++) {
      char group_name[80] = {0};
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
      group_new = H5Gopen   (file_id_,group_full.c_str());
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
  if (is_group_open_) {
    
    herr_t retval = H5Gclose(group_id_);

    ASSERT2("FileHdf5::group_close", "Return value %d closing group %s",
	    retval,group_name_.c_str(), (retval >= 0));

  }
  is_group_open_ = false;
}

//----------------------------------------------------------------------

void FileHdf5::group_read_meta
  ( void * buffer, std::string name,  scalar_type * type,
    int * nx, int * ny, int * nz) throw()
{
  // error check file open

  std::string file_name = path_ + "/" + name_;

  ASSERT1("FileHdf5::group_read_meta",
	 "Trying to read attribute from the unopened file %s",
	  file_name.c_str(),
	  is_file_open_);

  // error check group open

  ASSERT1("FileHdf5::group_read_meta",
	  "Trying to read attribute from unopened group %s",
	  group_name_.c_str(),
	  is_group_open_);

  hid_t meta_id = H5Aopen_name(group_id_, name.c_str());

  ASSERT3("FileHdf5::group_read_meta",
	  "H5Aopen_name() returned %d when opening attribute %s in file %s",
	   meta_id, name.c_str(),file_name.c_str(),
	   (meta_id >= 0));

  // Get attribute size

  hid_t meta_space_id = get_attr_space_ (meta_id,name);

  // set output extents

  get_output_extents_(meta_space_id,nx,ny,nz);

  // set output parameters

  scalar_type scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

  if (type) (*type) = scalar_type;

  // Read the attribute

  int retval = H5Aread
    (meta_id, scalar_to_hdf5_(scalar_type), buffer);

  // error check H5Aread

  ASSERT1("FileHdf5::group_read_meta","H5Aread returned %d",retval,(retval>=0));
}

//----------------------------------------------------------------------

void FileHdf5::set_compress (int level) throw ()
{
  compress_level_ = level; 
  if (compress_level_ != 0) {
    WARNING("FileHdf5::set_compress",
	    "Hard-coded for 2D data with chunk size [10,10]");
    int rank = 2;
    hsize_t chunk_size[MAX_DATA_RANK];
    chunk_size[0]=10;
    chunk_size[1]=10;
    H5Pset_chunk(data_prop_,rank,chunk_size);
    H5Pset_deflate(data_prop_,compress_level_);
  }
}

//======================================================================

void FileHdf5::write_meta_
( hid_t type_id,
  const void * buffer, std::string name, scalar_type type,
  int nx, int ny, int nz) throw()
{
  // error check file open


  std::string file_name = path_ + "/" + name_;

  ASSERT1("FileHdf5::write_meta_",
	  "Trying to write metadata to the unopened file %s",
	  file_name.c_str(), is_file_open_);

  // error check group or dataset open

  if (type_id == group_id_) {
    ASSERT1("FileHdf5::write_meta_",
	    "Trying to write attribute to unopened group %s",
	    group_name_.c_str(),is_group_open_);
  }
  if (type_id == data_id_) {
    ASSERT1("FileHdf5::write_meta",
	    "Trying to read attribute from unopened dataset %s",
	    data_name_.c_str(),is_data_open_);
  }

  // Determine the attribute rank

  // Create data space

  hid_t meta_space_id = create_data_space_ (nx,ny,nz,nx,ny,nz);

  // Create the attribute

  DEBUG1("Calling H5Acreate(%s)",name.c_str());
  hid_t meta_id = H5Acreate ( type_id,
			      name.c_str(),
			      scalar_to_hdf5_(type),
			      meta_space_id,
			      H5P_DEFAULT);

  // error check H5Acreate

  ASSERT2("FileHdf5::write_meta_","H5Acreate(%s) returned %d",
	  name.c_str(),meta_id,(meta_id>=0));

  // Write the attribute 

  DEBUG1("buffer = %p",buffer);
  H5Awrite (meta_id, scalar_to_hdf5_(type), buffer);

  // Close the attribute dataspace

  close_space_(meta_space_id);

  // Close the attribute

  int retval = H5Aclose(meta_id);

  ASSERT1("FileHdf5::write_meta_",
	  "H5Aclose() returned %d",retval,(retval >= 0));

}

//----------------------------------------------------------------------

int FileHdf5::scalar_to_hdf5_ (scalar_type type) const throw()
{
  // (*) NATIVE    -   FLOAT DOUBLE LDOUBLE
  // ( ) IEEE      -   F32BE F64BE     -
  // ( ) STD     B16BE B32BE B64BE     -
  // Types: http://www.hdfgroup.org/HDF5/Tutor/datatypes.html#native-types
  // char          H5T_NATIVE_CHAR   H5T_STD_I8BE   or H5T_STD_I8LE
  // float         H5T_NATIVE_FLOAT  H5T_IEEE_F32BE or H5T_IEEE_F32LE  
  // double        H5T_NATIVE_DOUBLE H5T_IEEE_F64BE or H5T_IEEE_F64LE  
  // unsigned char H5T_NATIVE_UCHAR  H5T_STD_U8BE   or H5T_STD_U8LE
  // int           H5T_NATIVE_INT    H5T_STD_I32BE  or H5T_STD_I32LE
  // short:        H5T_NATIVE_SHORT  H5T_STD_I16BE  or H5T_STD_I16LE
  // long:         H5T_NATIVE_LONG   H5T_STD_I32BE  or H5T_STD_I32LE
  //               H5T_STD_I64BE or  H5T_STD_I64LE
  // long long:    H5T_NATIVE_LLONG  H5T_STD_I64BE  or H5T_STD_I64LE

  hid_t hdf5_type;

  switch (type) {
  case scalar_type_unknown:
    ERROR("FileHdf5::scalar_to_hdf5_",
	  "scalar_type_unknown not implemented");
    hdf5_type = 0;
    break;
  case scalar_type_float:
    hdf5_type = H5T_NATIVE_FLOAT;
    break;
  case scalar_type_double:
    hdf5_type = H5T_NATIVE_DOUBLE;
    break;
  case scalar_type_long_double:
    ERROR("FileHdf5::scalar_to_hdf5_","long double not supported");
    hdf5_type = 0;
    break;
  case scalar_type_char:
    hdf5_type = H5T_NATIVE_CHAR;
    break;
  case scalar_type_int:
    hdf5_type = H5T_NATIVE_INT;
    break;
  case scalar_type_long:
    // use H5T_NATIVE_INT for long if same size: so inverse exists
    hdf5_type = (sizeof(long) == sizeof(int)) 
      ? H5T_NATIVE_INT : H5T_NATIVE_LONG;
    break;
  case scalar_type_long_long:
    hdf5_type = (sizeof(long long) == sizeof(long)) 
      ? H5T_NATIVE_LONG : H5T_NATIVE_LLONG;
    break;
  default:
    ERROR1("FileHdf5::scalar_to_hdf5_", "unsupported type %d", type);
    hdf5_type = 0;
    break;
  }
  return hdf5_type;
}

//----------------------------------------------------------------------

scalar_type FileHdf5::hdf5_to_scalar_ (int hdf5_type) const throw()
{

  H5T_class_t hdf5_class = H5Tget_class(hdf5_type);
  size_t      hdf5_size  = H5Tget_size (hdf5_type);

  scalar_type type = scalar_type_unknown;
 
  if (hdf5_class == H5T_INTEGER) {

    if (hdf5_size == sizeof(char)) {
      type = scalar_type_char;
    } else if (hdf5_size == sizeof(int)) {
      type = scalar_type_int;
    } else if (hdf5_size == sizeof(long)) {
      type = scalar_type_long;
    } else if (hdf5_size == sizeof(long long)) {
      type = scalar_type_long_long;
    } else ASSERT("","",0);

  } else if (hdf5_class == H5T_FLOAT) {

    if (hdf5_size == sizeof(float)) {
      type = scalar_type_float;
    } else if (hdf5_size == sizeof(double)) {
      type = scalar_type_double;
    } else ASSERT("","",0);

  } else {

    ERROR2("FileHdf5::hdf5_to_scalar_",
	  "Unknown type of class %d and size %d",
	  hdf5_class, int(hdf5_size));
  }

  return type;

}

//----------------------------------------------------------------------

std::string FileHdf5::relative_to_absolute_  
(
 std::string path_relative, 
 std::string path_absolute
 ) const throw()
{
  // First add trailing "/" if needed

  if (path_absolute[path_absolute.size()-1] != '/') {
    path_absolute = path_absolute + "/";
  }
    
  // Return "relative" path if it's already absolute

  if (path_relative[0] == '/') return path_relative;

  // Start with existing absolute path, and construct new absolute path
  // given relative path

  std::string path_dir;
  
  size_t p_left_slash=path_relative.find("/");

  while (p_left_slash != std::string::npos) {

    path_dir      = path_relative.substr(0,p_left_slash);
    path_relative = path_relative.substr(p_left_slash+1,std::string::npos);

    if (path_dir == "..") {
      int len= path_absolute.size();
      int p_right_slash = path_absolute.rfind ("/",len-2);
      path_absolute = path_absolute.substr (0,p_right_slash+1);
    } else {
      path_absolute = path_absolute + path_dir + "/";
    }

    p_left_slash=path_relative.find("/");

  }
  path_dir = path_relative;
  if (path_dir == "..") {
    int len= path_absolute.size();
    int p_right_slash = path_absolute.rfind ("/",len-2);
    path_absolute = path_absolute.substr (0,p_right_slash+1);
  } else {
    path_absolute = path_absolute + path_dir + "/";
  }
  
  // remove trailing "/" if needed
  if (path_absolute.size() > 1) {
    path_absolute = path_absolute.substr(0,path_absolute.size()-1);
  }
  return path_absolute;
}

//----------------------------------------------------------------------

void FileHdf5::get_output_extents_
( hid_t data_space_id, int * nx, int * ny, int * nz) throw ()
{

   hsize_t data_size[MAX_DATA_RANK];
   int rank = H5Sget_simple_extent_dims(data_space_id,data_size,0);

  // error check rank

  ASSERT1("FileHdf5::get_output_extents_","rank %d is out of range",
	  rank, (1 <= rank && rank <= MAX_DATA_RANK));

   // Get data size: NOTE REVERSED AXES

  if (rank == 1) {
    if (nx) (*nx) = data_size[0];
    if (ny) (*ny) = 1;
    if (nz) (*nz) = 1;
  }
  if (rank == 2) {
    if (nx) (*nx) = data_size[1];
    if (ny) (*ny) = data_size[0];
    if (nz) (*nz) = 1;
  }
  if (rank == 3) {
    if (nx) (*nx) = data_size[2];
    if (ny) (*ny) = data_size[1];
    if (nz) (*nz) = data_size[0];
  }
}

//----------------------------------------------------------------------

hid_t FileHdf5::create_space_(int nxd,int nyd,int nzd,
			      int nx,int ny,int nz) throw ()
{

  hsize_t space_dims[MAX_DATA_RANK];
  hsize_t space_size[MAX_DATA_RANK];

  hsize_t rank = MAX_DATA_RANK;

  if (nzd == 0 || nzd == 1) -- rank;
  if (nyd == 0 || nyd == 1) -- rank;


  // error check rank

  ASSERT1("FileHdf5::create_space_","rank %d is out of range",
	  rank, (1 <= rank && rank <= MAX_DATA_RANK));

  // Define the space: NOTE REVERSED AXES

  bool need_dims = true;

  if (rank == 1) {

    space_dims[0] = nxd;

    space_size[0] = nx ? nx : nxd;  

    need_dims = (space_dims[0] != space_size[0]);

  }
  if (rank == 2) {

    space_dims[0] = nyd;
    space_dims[1] = nxd;

    space_size[0] = ny ? ny : nyd;  
    space_size[1] = nx ? nx : nxd;  

    need_dims = (space_dims[0] != space_size[0] ||
		 space_dims[1] != space_size[1]);

  }
  if (rank == 3) {

    space_dims[0] = nzd;
    space_dims[1] = nyd;
    space_dims[2] = nxd;

    space_size[0] = nz ? nz : nzd;
    space_size[1] = ny ? ny : nyd;
    space_size[2] = nx ? nx : nxd;

    need_dims = (space_dims[0] != space_size[0] ||
		 space_dims[1] != space_size[1] ||
		 space_dims[2] != space_size[2]);

  }

  hid_t space_id;

  if (need_dims) {

    space_id = H5Screate_simple (rank, space_dims, 0);

    hsize_t start[3];
    hsize_t count[3];

    for (size_t i=0; i<rank; i++) {
      start[i] = (space_dims[i]-space_size[i])/2;
      count[i] = space_size[i];
    }

    H5Sselect_hyperslab (space_id,H5S_SELECT_SET,start,0,count,0);
  }
  else space_id = H5Screate_simple (rank, space_dims, 0); 

  // error check H5Screate_simple

  ASSERT1("FileHdf5::create_space_",
	  "h5Screate_simple returned %d",
	  space_id,
	  (space_id>=0));

  return space_id;
}

//----------------------------------------------------------------------

void FileHdf5::close_space_ (hid_t space_id) throw()
{
  // Close space

  int retval = 0;

  if (space_id != H5S_ALL) retval = H5Sclose (space_id);

  // Error check H5Sclose 

  ASSERT1("FileHdf5::close_space_", "Return value %d", retval,(retval >= 0));
}

//----------------------------------------------------------------------

hid_t FileHdf5::open_dataset_ (hid_t group, std::string name) throw()
{
  hid_t dataset_id = H5Dopen( group, name.c_str());

  // error check H5Dopen

  ASSERT3("FileHdf5::open_dataset_", 
	  "H5Dopen() returned %d opening %s in file %s",
	  dataset_id,name.c_str(),(path_ + "/" + name_).c_str(), 
	  data_id_ >= 0);

  return dataset_id;

}
//----------------------------------------------------------------------

void FileHdf5::close_dataset_ () throw()
{
  int retval = H5Dclose (data_id_);

  // error check H5Dclose

  ASSERT2("FileHdf5::close_dataset_", "Return value %d closing dataset %s",
	  retval,data_name_.c_str(), (retval >= 0));
}

//----------------------------------------------------------------------

hid_t FileHdf5::get_data_space_(hid_t data_id, std::string name) throw ()
{

  hid_t data_space_id = H5Dget_space (data_id);

  // error check rank

  int rank = H5Sget_simple_extent_ndims(data_space_id);

  ASSERT3("FileHdf5::get_data_space_", 
	  "Dataset %s in file %s has unsupported rank %d",
	  name.c_str(),(path_ + "/" + name_).c_str(),rank,
	  (1 <= rank && rank <= MAX_DATA_RANK));

  return data_space_id;

}
//----------------------------------------------------------------------

hid_t FileHdf5::get_attr_space_(hid_t attr_id, std::string name) throw ()
{

  hid_t data_space_id = H5Aget_space (attr_id);

  // error check rank

  int rank = H5Sget_simple_extent_ndims(data_space_id);

  ASSERT3("FileHdf5::get_attr_space_",
	  "Attribute %s in file %s has unsupported rank %d",
	  name.c_str(),(path_ + "/" + name_).c_str(),rank,
	  (1 <= rank && rank <= MAX_DATA_RANK));

  return data_space_id;

}
