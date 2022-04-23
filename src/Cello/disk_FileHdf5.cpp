// See LICENSE_CELLO file for license and copyright information

/// @file      disk_FileHdf5.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:11:36 PST 2008
/// @brief     Implementation of the FileHdf5 class

#include "cello.hpp"

#include "disk.hpp"

// #define TRACE_DISK

#define MAX_DATA_RANK 4
#define MAX_ATTR_RANK 4

//----------------------------------------------------------------------

std::map<const std::string,FileHdf5 *> FileHdf5::file_list;

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
    data_type_(type_unknown),
    data_rank_(0),
    data_prop_(H5P_DEFAULT),
    is_data_open_(false),
    compress_level_(0)
{
  for (int i=0; i<MAX_DATA_RANK; i++) {
    data_dims_[i] = 0;
  }

  // data_prop_ = H5P_DEFAULT;
  // group_prop_ = H5Pcreate (H5P_GROUP_CREATE);

  data_prop_  = H5Pcreate (H5P_DATASET_CREATE);
#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Pcreate(%d)\n",CkMyPe(),file_id_, __LINE__,data_prop_);
  fflush(stdout);
#endif  
  group_prop_ = H5P_DEFAULT;
}

//----------------------------------------------------------------------

FileHdf5::~FileHdf5() throw()
{
#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Pclose(%d)\n",CkMyPe(),file_id_, __LINE__,data_prop_);
  fflush(stdout);
#endif  
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
#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Fopen(%s)\n",CkMyPe(),file_id_, __LINE__,file_name.c_str());
  fflush(stdout);
#endif  

  // error check file opened

  ASSERT2("FileHdf5::file_open", "Return value %ld opening file %s",
	 file_id_,file_name.c_str(), file_id_ >= 0);

  // update file state
  is_file_open_ = true;

}

//----------------------------------------------------------------------

int FileHdf5::data_size (int * m4_int) throw()
{
  hsize_t m4[4] = {0};
  hsize_t n4[4] = {0};
  int rank = H5Sget_simple_extent_dims(data_space_id_,n4,m4);
  m4_int[0] = int(m4[0]);
  m4_int[1] = int(m4[1]);
  m4_int[2] = int(m4[2]);
  m4_int[3] = int(m4[3]);
  return rank;
}

//----------------------------------------------------------------------

void FileHdf5::file_create () throw()
{

  // create file

  std::string file_name = path_ + "/" + name_;

  file_id_ = H5Fcreate(file_name.c_str(),
		       H5F_ACC_TRUNC,
		       H5P_DEFAULT,
		       H5P_DEFAULT);
#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Fcreate(%d)\n",CkMyPe(),file_id_, __LINE__,file_id_);
  fflush(stdout);
#endif  

  // error check file created

  TRACE1("File = %s",file_name.c_str());
  ASSERT2("FileHdf5::file_create",  "Return value %ld opening file %s",
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
  data_close();

  // Close the file

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Fclose(%s)\n",CkMyPe(),file_id_, __LINE__,file_name.c_str());
  fflush(stdout);
#endif  
  int retval = H5Fclose (file_id_);

  // error check H5Fclose

  ASSERT2("FileHdf5::file_close", "Return value %d closing file %s",
	  retval,file_name.c_str(), (retval >= 0));

  // update file state

  is_file_open_ = false;
}

//----------------------------------------------------------------------

void FileHdf5::data_open
( std::string name,  int * type,
  int * m1, int * m2, int * m3, int * m4) throw()
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

  get_extents_(data_space_id_,m1,m2,m3,m4);

  // find subset if selecting subset

  // Initialize name and type

  data_name_ = name;
  data_type_ = hdf5_to_scalar_(H5Dget_type (data_id_));
  //  WARNING("FileHdf5::data_open","Set memspace");

  // set output parameters

  if (type) (*type) = data_type_;

}

//----------------------------------------------------------------------

void FileHdf5::data_slice
( int m1, int m2, int m3, int m4,
  int n1, int n2, int n3, int n4,
  int o1, int o2, int o3, int o4) throw()
{
  data_space_id_ = space_slice_ (data_space_id_,
				 m1, m2, m3,m4,
				 n1,n2,n3,n4,
				 o1,o2,o3,o4);
}

//----------------------------------------------------------------------

void FileHdf5::mem_create
( int mx, int my, int mz,
  int nx, int ny, int nz,
  int gx, int gy, int gz )
{
  mem_space_id_ = space_create_ (mx,my,mz,1, nx,ny,nz,1, gx,gy,gz,0);
}


//----------------------------------------------------------------------

void FileHdf5::data_create
( std::string name,  int type,
  int m1, int m2, int m3, int m4,
  int n1, int n2, int n3, int n4,
  int o1, int o2, int o3, int o4) throw()
{

  if (n1==0) n1=m1;
  if (n2==0) n2=m2;
  if (n3==0) n3=m3;
  if (n4==0) n4=m4;

  // Initialize data attributes
  data_name_ = name;
  data_type_ = type;

  // Determine the group

  hid_t group = (is_group_open_) ? group_id_ : file_id_;

  // Create dataspace

  data_space_id_ = space_create_ (m1,m2,m3,m4,
				  n1,n2,n3,n4,
				  o1,o2,o3,o4);

  // Create the new dataset

  data_id_ = H5Dcreate( group,
			name.c_str(),
			scalar_to_hdf5_(type),
			data_space_id_,
			H5P_DEFAULT,
			data_prop_,
			H5P_DEFAULT);
#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Dcreate(%d)\n",CkMyPe(),file_id_, __LINE__,data_id_);
  fflush(stdout);
#endif  

  // error check H5Dcreate

  ASSERT2("FileHdf5::data_create", "Return value %ld creating dataset %s",
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

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Dread(%d)\n",CkMyPe(),file_id_, __LINE__,data_id_);
  fflush(stdout);
#endif  
  int retval = 
    H5Dread (data_id_,
	     scalar_to_hdf5_(data_type_),
	     mem_space_id_,
	     data_space_id_,
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

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Dwrite(%d)\n",CkMyPe(),file_id_, __LINE__,data_id_);
  fflush(stdout);
#endif  
  int retval = 
    H5Dwrite (data_id_,
	      scalar_to_hdf5_(data_type_),
	      mem_space_id_,
	      data_space_id_,
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

    space_close_(data_space_id_);
    //     space_close_(mem_space_id_);

    // close the dataset

    close_dataset_ ();

    // update dataset state

    is_data_open_ = false;
  }
}

//----------------------------------------------------------------------

void FileHdf5::file_read_scalar
( void * buffer, std::string name,  int * type) throw()
{
  std::string file_name = path_ + "/" + name_;

  // error check file open

  ASSERT1("FileHdf5::file_read_meta",
	  "Trying to read metadata from the unopened file %s",
	  file_name.c_str(), is_file_open_);

  hid_t meta_id = H5Aopen_name(file_id_, name.c_str());

  // error check H5Aopen_name

  ASSERT3("FileHdf5::file_read_meta",
	  "H5Aopen_name() returned %ld opening %s in file %s",
	  meta_id, name.c_str(),file_name.c_str(),
	  (meta_id >= 0));

  // set output type

  int scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

  if (type) (*type) = scalar_type;

  // Read the attribute

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Aread()\n",CkMyPe(),file_id_, __LINE__);
  fflush(stdout);
#endif  
  int retval = 
    H5Aread(meta_id, scalar_to_hdf5_(scalar_type), buffer);

  // error check H5Aread

  ASSERT1("FileHdf5::file_read_meta_","H5Aread() returned %d",
	  retval,(retval>=0));
}


//----------------------------------------------------------------------

void FileHdf5::file_read_meta
  ( void * buffer, std::string name,  int * type,
    int * n1, int * n2, int * n3, int * n4) throw()
{


  std::string file_name = path_ + "/" + name_;

  // error check file open

  ASSERT1("FileHdf5::file_read_meta",
	  "Trying to read metadata from the unopened file %s",
	  file_name.c_str(), is_file_open_);

  hid_t meta_id = H5Aopen_name(file_id_, name.c_str());

  // error check H5Aopen_name

  ASSERT3("FileHdf5::file_read_meta",
	  "H5Aopen_name() returned %ld opening %s in file %s",
	  meta_id, name.c_str(),file_name.c_str(),
	  (meta_id >= 0));

  // get dataspace
  hid_t meta_space_id = get_attr_space_(meta_id,name);

  // set output extents

  get_extents_(meta_space_id,n1,n2,n3,n4);

  // set output type

  int scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

  if (type) (*type) = scalar_type;

  // Read the attribute

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Aread()\n",CkMyPe(),file_id_, __LINE__);
  fflush(stdout);
#endif  
  int retval = 
    H5Aread(meta_id, scalar_to_hdf5_(scalar_type), buffer);

  // error check H5Aread

  ASSERT1("FileHdf5::file_read_meta_","H5Aread() returned %d",
	  retval,(retval>=0));

}

//----------------------------------------------------------------------

void FileHdf5::data_read_meta
  ( void * buffer, std::string name,  int * type,
    int * n1, int * n2, int * n3, int * n4) throw()
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
	  "H5Aopen_name() returned %ld when opening attribute %s in file %s",
	   meta_id, name.c_str(),file_name.c_str(),
	   (meta_id >= 0));

  // Get attribute size

  hid_t meta_space_id = get_attr_space_ (meta_id,name);

  // set output extents

  get_extents_(meta_space_id,n1,n2,n3,n4);

  // set output parameters

  int scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

  if (type) (*type) = scalar_type;

  // Read the attribute

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Aread()\n",CkMyPe(),file_id_, __LINE__);
  fflush(stdout);
#endif  
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

  group_id_ = H5Gopen(file_id_, group_name_.c_str(),H5P_DEFAULT);
#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Gopen(%d)\n",CkMyPe(),file_id_, __LINE__,group_id_);
  fflush(stdout);
#endif  

  // error check H5Gopen()

  ASSERT2("FileHdf5::group_open()", "H5Gopen(%s) returned %ld", 
	  group_name_.c_str(),group_id_,group_id_>=0);

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

  group_id_ = H5Gopen(file_id_,group_full.c_str(),H5P_DEFAULT);
#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Gopen(%d)\n",CkMyPe(),file_id_, __LINE__,group_id_);
  fflush(stdout);
#endif  
  ASSERT2("FileHdf5::group_open()", "H5Gopen(%s) returned %ld", 
	  group_full.c_str(),group_id_,group_id_>=0);

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
    H5Gget_info_by_name (file_id_, group_full.c_str(),
			 &group_info, H5P_DEFAULT);

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
      group_new = H5Gopen   (file_id_,group_full.c_str(), H5P_DEFAULT);
#ifdef TRACE_DISK  
      CkPrintf ("%d %Ld :%d TRACE_DISK H5Gopen(%d)\n",CkMyPe(),file_id_, __LINE__,group_new);
  fflush(stdout);
#endif  
    } else {
      group_new = H5Gcreate (file_id_,group_full.c_str(),
			     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#ifdef TRACE_DISK  
      CkPrintf ("%d %Ld :%d TRACE_DISK H5Gcreate(%d)\n",CkMyPe(),file_id_, __LINE__,group_new);
  fflush(stdout);
#endif  
    }

    // Close parent group

#ifdef TRACE_DISK  
    CkPrintf ("%d %Ld :%d TRACE_DISK H5Gclose(%d)\n",CkMyPe(),file_id_, __LINE__,group_id_);
  fflush(stdout);
#endif  
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

#ifdef TRACE_DISK  
    CkPrintf ("%d %Ld :%d TRACE_DISK H5Gclose(%d)\n",CkMyPe(),file_id_, __LINE__,group_id_);
  fflush(stdout);
#endif  

    herr_t retval = H5Gclose(group_id_);

    ASSERT2("FileHdf5::group_close", "Return value %d closing group %s",
	    retval,group_name_.c_str(), (retval >= 0));

  }
  is_group_open_ = false;
}

//----------------------------------------------------------------------

void FileHdf5::group_read_meta
  ( void * buffer, std::string name,  int * type,
    int * n1, int * n2, int * n3, int * n4) throw()
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
	  "H5Aopen_name() returned %ld when opening attribute %s in file %s",
	   meta_id, name.c_str(),file_name.c_str(),
	   (meta_id >= 0));

  // Get attribute size

  hid_t meta_space_id = get_attr_space_ (meta_id,name);

  // set output extents

  get_extents_(meta_space_id,n1,n2,n3,n4);

  // set output parameters

  int scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

  if (type) (*type) = scalar_type;

  // Read the attribute

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Aread()\n",CkMyPe(),file_id_, __LINE__);
  fflush(stdout);
#endif  
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
  const void * buffer, std::string name, int type,
  int n1, int n2, int n3, int n4) throw()
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

  hid_t meta_space_id = space_create_ (n1,n2,n3,n4,n1,n2,n3,n4,0,0,0,0);

  // Create the attribute

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Acreate(%s)\n",CkMyPe(),file_id_, __LINE__,name.c_str());
  fflush(stdout);
#endif  
  hid_t meta_id = H5Acreate ( type_id,
			      name.c_str(),
			      scalar_to_hdf5_(type),
			      meta_space_id,
			      H5P_DEFAULT,
			      H5P_DEFAULT);

  // error check H5Acreate

  ASSERT2("FileHdf5::write_meta_","H5Acreate(%s) returned %ld",
	  name.c_str(),meta_id,(meta_id>=0));

  // Write the attribute 

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Awrite()\n",CkMyPe(),file_id_, __LINE__);
  fflush(stdout);
#endif  
  H5Awrite (meta_id, scalar_to_hdf5_(type), buffer);

  // Close the attribute dataspace

  space_close_(meta_space_id);

  // Close the attribute

#ifdef TRACE_DISK  
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Aclose()\n",CkMyPe(),file_id_, __LINE__);
  fflush(stdout);
#endif  
  int retval = H5Aclose(meta_id);

  ASSERT1("FileHdf5::write_meta_",
	  "H5Aclose() returned %d",retval,(retval >= 0));

}

//----------------------------------------------------------------------

hid_t FileHdf5::scalar_to_hdf5_ (int type) const throw()
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

  // Whether to use native encodings.  Note, setting to false doesn't
  // work--written values are garbage when viewed using h5dump
  const bool native = true;

  // Whether to use big-endian or little-endian.  Ignored if native
  // is true
  bool be = true; 

  switch (type) {
  case type_unknown:
  case type_extended80:
  case type_extended96:
    ERROR1("FileHdf5::scalar_to_hdf5_",
	   "type_%s not implemented",cello::type_name[type]);
    hdf5_type = 0;
    break;
  case type_single:
    hdf5_type = native ? 
      H5T_NATIVE_FLOAT : (be ? H5T_IEEE_F32BE : H5T_IEEE_F32LE);
    break;
  case type_double:
    hdf5_type = native ?
      H5T_NATIVE_DOUBLE : (be ? H5T_IEEE_F64BE : H5T_IEEE_F64LE);
    break;
  case type_quadruple:
    if (native) {
      hdf5_type = H5T_NATIVE_LDOUBLE;
    } else {
      ERROR("FileHdf5::scalar_to_hdf5_",
	    "type_quadruple not implemented unless native=true");
    }
    break;
  case type_char:
    hdf5_type = H5T_NATIVE_CHAR;
    break;
  case type_short:
    hdf5_type = native ? 
      H5T_NATIVE_SHORT : (be ? H5T_STD_I16BE : H5T_STD_I16LE);
    break;
  case type_int:
    hdf5_type = native ? 
      H5T_NATIVE_INT : (be ? H5T_STD_I32BE : H5T_STD_I32LE);
    break;
  case type_long_long:
    hdf5_type = native ? 
      H5T_NATIVE_LLONG : (be ? H5T_STD_I64BE : H5T_STD_I64LE);
    break;
  default:
    ERROR1("FileHdf5::scalar_to_hdf5_", "unsupported type %d", type);
    hdf5_type = 0;
    break;
  }
  return hdf5_type;
}

//----------------------------------------------------------------------

int FileHdf5::hdf5_to_scalar_ (hid_t hdf5_type) const throw()
{

  H5T_class_t hdf5_class = H5Tget_class(hdf5_type);
  size_t      hdf5_size  = H5Tget_size (hdf5_type);

  int type = type_unknown;
 
  if (hdf5_class == H5T_INTEGER) {

    if (hdf5_size == sizeof(char)) {
      type = type_char;
    } else if (hdf5_size == sizeof(int)) {
      type = type_int;
    } else if (hdf5_size == sizeof(long long)) {
      type = type_long_long;
    } else ASSERT("","",0);

  } else if (hdf5_class == H5T_FLOAT) {

    if (hdf5_size == sizeof(float)) {
      type = type_float;
    } else if (hdf5_size == sizeof(double)) {
      type = type_double;
    } else if (hdf5_size == sizeof(long double)) {
      type = type_quadruple;
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

void FileHdf5::get_extents_
( hid_t data_space_id, int * n1, int * n2, int * n3, int * n4) throw ()
{

   hsize_t data_size[MAX_DATA_RANK];
   int rank = H5Sget_simple_extent_dims(data_space_id,data_size,0);

  // error check rank

  ASSERT1("FileHdf5::get_extents_","rank %d is out of range",
	  rank, (1 <= rank && rank <= MAX_DATA_RANK));

   // Get data size: NOTE REVERSED AXES

  if (rank == 1) {
    if (n1) (*n1) = data_size[0];
    if (n2) (*n2) = 1;
    if (n3) (*n3) = 1;
    if (n4) (*n4) = 1;
  }
  if (rank == 2) {
    if (n1) (*n1) = data_size[0];
    if (n2) (*n2) = data_size[1];
    if (n3) (*n3) = 1;
    if (n4) (*n4) = 1;
  }
  if (rank == 3) {
    if (n1) (*n1) = data_size[0];
    if (n2) (*n2) = data_size[1];
    if (n3) (*n3) = data_size[2];
    if (n4) (*n4) = 1;
  }
  if (rank == 4) {
    if (n1) (*n1) = data_size[0];
    if (n2) (*n2) = data_size[1];
    if (n3) (*n3) = data_size[2];
    if (n4) (*n4) = data_size[3];
  }
}

//----------------------------------------------------------------------

hid_t FileHdf5::space_create_(int m1, int m2, int m3, int m4,
			      int n1, int n2, int n3, int n4,
			      int o1, int o2, int o3, int o4) throw ()
{
  hsize_t dims[4] = {1};
  hsize_t count[4] = {1};
  hsize_t start[4] = {0};

  hsize_t rank = 4;

  if (m4 == 0 || m4 == 1) -- rank;
  if (m3 == 0 || m3 == 1) -- rank;
  if (m2 == 0 || m2 == 1) -- rank;

  // error check rank

  ASSERT1("FileHdf5::space_create_","rank %llu is out of range",
	  rank, (1 <= rank && rank <= 4));

  // Define the space: NOTE REVERSED AXES

  bool need_hyper = true;

  if (rank == 1) {

    need_hyper = ((n1 != m1) || (o1 != 0));

    dims[0] = m1;

    count[0] = n1 ? n1 : m1;

    start[0] = o1;

  } else if (rank == 2) {

    need_hyper = ((n1 != m1) || (o1 != 0) ||
		  (n2 != m2) || (o2 != 0));
    
    dims[0] = m1;
    dims[1] = m2;

    count[0] = n1 ? n1 : m1;
    count[1] = n2 ? n2 : m2;  

    start[0] = o1;
    start[1] = o2;

  } else if (rank == 3) {

    need_hyper = ((n1 != m1) || (o1 != 0) ||
		  (n2 != m2) || (o2 != 0) ||
		  (n3 != m3) || (o3 != 0));

    dims[0] = m1;
    dims[1] = m2;
    dims[2] = m3;

    count[0] = n1 ? n1 : m1;
    count[1] = n2 ? n2 : m2;
    count[2] = n3 ? n3 : m3;

    start[0] = o1;
    start[1] = o2;
    start[2] = o3;

  } else if (rank == 4) {

    need_hyper = ((n1 != m1) || (o1 != 0) ||
		  (n2 != m2) || (o2 != 0) ||
		  (n3 != m3) || (o3 != 0) ||
		  (n4 != m4) || (o4 != 0));

    dims[0] = m1;
    dims[1] = m2;
    dims[2] = m3;
    dims[3] = m4;

    count[0] = n1 ? n1 : m1;
    count[1] = n2 ? n2 : m2;
    count[2] = n3 ? n3 : m3;
    count[3] = n4 ? n4 : m4;

    start[0] = o1;
    start[1] = o2;
    start[2] = o3;
    start[3] = o4;

  } else {
    ERROR1 ("FileHdf5::space_create_()",
	    "Rank %llu out of bounds 1 <= rank <= 4",rank);
  }

  hid_t space_id = H5Screate_simple (rank, dims, 0);

  if (need_hyper) {
    H5Sselect_hyperslab (space_id,H5S_SELECT_SET,start,0,count,0);
  }

  // error check H5Screate_simple

  ASSERT1("FileHdf5::space_create_",
	  "h5Screate_simple returned %ld",
	  space_id,
	  (space_id>=0));

  return space_id;
}

//----------------------------------------------------------------------

hid_t FileHdf5::space_slice_(hid_t space_id,
			     int m1, int m2, int m3, int m4,
			     int n1, int n2, int n3, int n4,
			     int o1, int o2, int o3, int o4) throw ()
{
  int rank = H5Sget_simple_extent_ndims(space_id);
  bool need_hyper = false;
  int M4[4] = {m1,m2,m3,m4};
  int N4[4] = {n1,n2,n3,n4};
  int O4[4] = {o1,o2,o3,o4};
  for (int i=0; i<rank; i++) {
    if (O4[i]!=0 || N4[i] != M4[i]) need_hyper = true;
  }
  if (! need_hyper) return space_id;

  hsize_t start[4]  = {hsize_t(o1),
		       hsize_t(o2),
		       hsize_t(o3),
		       hsize_t(o4)};
  hsize_t count[4]  = {hsize_t(n1),
		       hsize_t(n2),
		       hsize_t(n3),
		       hsize_t(n4)};

  H5Sselect_hyperslab (space_id,H5S_SELECT_SET,start,0,count,0);

  return space_id;
}

//----------------------------------------------------------------------

void FileHdf5::space_close_ (hid_t space_id) throw()
{
  // Close space

  int retval = 0;

  if (space_id != H5S_ALL) retval = H5Sclose (space_id);

  // Error check H5Sclose 

  ASSERT1("FileHdf5::space_close_", "Return value %d", retval,(retval >= 0));
}

//----------------------------------------------------------------------

hid_t FileHdf5::open_dataset_ (hid_t group, std::string name) throw()
{
  
#ifdef TRACE_DISK
  CkPrintf ("%d %Ld :%d TRACE_DISK H5Dopen(%s)\n",CkMyPe(),file_id_, __LINE__,name.c_str());
  fflush(stdout);
#endif  
  hid_t dataset_id = H5Dopen( group, name.c_str(), H5P_DEFAULT);

  // error check H5Dopen

  ASSERT3("FileHdf5::open_dataset_", 
	  "H5Dopen() returned %ld opening %s in file %s",
	  dataset_id,name.c_str(),(path_ + "/" + name_).c_str(), 
	  dataset_id >= 0);

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
