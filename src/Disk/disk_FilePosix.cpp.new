// // See LICENSE_CELLO file for license and copyright information

// /// @file      disk_FilePosix.cpp
// /// @author    James Bordner (jobordner@ucsd.edu)
// /// @date      2011-09-15
// /// @brief     Implementation of the FilePosix class

// #include "cello.hpp"

// #include "disk.hpp"
 
// //----------------------------------------------------------------------

// FilePosix::FilePosix (std::string path, std::string name) throw()
//   : File(path,name),
//     fp_(0),
//     is_file_open_(false),
//     data_type_(scalar_type_unknown),
//     data_rank_(0),
//     is_data_open_(false)
// {
//   const int rank_max = 5;
//   for (int i=0; i<rank_max; i++) {
//     data_size_[i] = 0;
//   }
// }

// //----------------------------------------------------------------------

// void FilePosix::file_open () throw()
// {
//   if (is_file_open_) {

//     WARNING1("FilePosix::file_open",
// 	    "Attempting to reopen an opened file %s",
// 	    name_.c_str());

//   } else {

//     std::string full_name = path_ + "/" + name_;

//     fp_ = fopen (full_name.c_str(), "rb");

//     // error check

//     if ( fp_ != NULL) {

//       is_file_open_ = true;

//     } else {

//       WARNING1("FilePosix::file_open",
// 	       "Error opening file %s",
// 	       full_name.c_str());

//     }
//   }
// }

// //----------------------------------------------------------------------

// void FilePosix::file_create () throw()
// {

//   std::string full_name = path_ + "/" + name_;

//   fp_ = fopen (full_name.c_str(), "wb");

//   // error check

//   if (fp_ != 0) {

//     is_file_open_ = true;

//   } else {

//     WARNING2("FilePosix::file_create",
// 	    "Return value %d opening file %s",
// 	    file_id_,full_name.c_str());

//   }
// }

// //----------------------------------------------------------------------

// void FilePosix::file_close () throw()
// {

//   int retval = fclose (fp_);

//   if (retval == 0) {
    
//     fp_ = 0;

//     is_file_open_ = false;

//   } else {
    
//     WARNING3("FilePosix::file_close",
// 	     "Return value %d closing file %s/%s",
// 	     retval,path_.c_str(),name_.c_str());
//   }
// }

// //----------------------------------------------------------------------

// void FilePosix::data_open
// ( std::string name,  enum scalar_type * type,
//   int * n0, int * n1, int * n2, int * n3, int * n4) throw()
// {

//  // Error checking

//   if (! is_file_open_) {

//     ERROR1("FilePosix::data_open",
// 	  "Trying to open dataset in unopened file %s",
// 	  (path_ + name_).c_str());
//   }

//   // Return 0: must use data_open_meta for type and dimensions if needed

//   if (type) (*type) = 0;

//   if (n0) (*n0) = 0;
//   if (n1) (*n1) = 0;
//   if (n2) (*n2) = 0;
//   if (n3) (*n3) = 0;
//   if (n4) (*n4) = 0;

// }

// //----------------------------------------------------------------------

// void FilePosix::data_create
// ( std::string name,  enum scalar_type type,
//   int n0, int n1, int n2, int n3, int n4) throw()
// {

//   // Initialize data attributes

//   data_name_ = name;
//   data_type_ = type;

//   data_rank_ = 5;

//   if (n4 == 0) -- data_rank_;
//   if (n3 == 0) -- data_rank_;
//   if (n2 == 0) -- data_rank_;
//   if (n1 == 0) -- data_rank_;

//   ASSERT ("FilePosix::data_set","d is out of range",
// 	  (1 <= data_rank_ && data_rank_ <= 5));

//   data_size_[0] = n0;
//   data_size_[1] = n1;
//   data_size_[2] = n2;
//   data_size_[3] = n3;
//   data_size_[4] = n4;


//   data_type_ = type;

//   is_data_open_ = true;


// }

// //----------------------------------------------------------------------

// void FilePosix::data_read
// ( void * buffer) throw()
// {

//   // error check file open

//   if (! is_file_open_) {

//     ERROR1("FilePosix::data_read",
// 	  "Trying to read from unopened file %s",
// 	  (path_ + name_).c_str());
//   }

//   // error check dataset open

//   if (! is_data_open_) {

//     ERROR1("FilePosix::data_read",
// 	  "Trying to read unopened dataset %s",
// 	  data_name_.c_str());
//   }

//   // read data

//   int bytes =  data_size_[0] 
//     *          data_size_[1]
//     *          data_size_[2]
//     *          data_size_[3]
//     *          data_size_[4];

//   bytes *= scalar_size_ (data_type_);

//   read(fp_,buffer, bytes);

// }

// //----------------------------------------------------------------------
// // @@@@@@@@@@@@@@@@@@@@@@@@

// void FilePosix::data_write
// ( const void * buffer) throw()
// {
//   // error check file open

//   if (! is_file_open_) {

//     ERROR1("FilePosix::data_write",
// 	  "Trying to write to unopened file %s",
// 	  (path_ + name_).c_str());
//   }

//   // error check dataset open

//   if (! is_data_open_) {

//     ERROR1("FilePosix::data_write",
// 	  "Trying to write unopened dataset %s",
// 	  data_name_.c_str());
//   }

//   // ...

// }

// //----------------------------------------------------------------------

// void FilePosix::data_close() throw()
// {
//   if (is_data_open_) {
//     fclose_ (fp_);
//     is_data_open_ = false;
//   }
// }

// //----------------------------------------------------------------------

// void FilePosix::file_read_meta
//   ( void * buffer, std::string name,  enum scalar_type * type,
//     int * n0, int * n1, int * n2, int * n3, int * n4) throw()
// {
//   if ( ! is_file_open_) {

//     ERROR1("FilePosix::file_read_meta",
// 	  "Trying to read metadata from the unopened file %s",
// 	  (path_ + name_).c_str());
//   }

//   // Get attribute size

//   hid_t meta_space_id = H5Aget_space (meta_id);

//   // error check rank

//   int rank = H5Sget_simple_extent_ndims(meta_space_id);

//   if (rank > 5) {
//     ERROR3("FilePosix::file_read_meta",
// 	  "Attribute %s in file %s has unsupported rank %d",
// 	  name.c_str(),name_.c_str(),rank);
//   }

//   // set output parameters

//   scalar_type scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

//   if (type) (*type) = scalar_type;

//   hsize_t n[5];
//   H5Sget_simple_extent_dims(meta_space_id,n,0);

//   if (n0) (*n0) = n[0];
//   if (n1) (*n1) = n[1];
//   if (n2) (*n2) = n[2];
//   if (n3) (*n3) = n[3];
//   if (n4) (*n4) = n[4];

//   // Read the attribute

//   H5Aread(meta_id, scalar_to_hdf5_(scalar_type), buffer);

// }

// //----------------------------------------------------------------------

// void FilePosix::file_write_meta
//   ( const void * buffer, std::string name, enum scalar_type type,
//     int n0, int n1, int n2, int n3, int n4) throw()
// {
//   if ( ! is_file_open_) {
//     ERROR1("FilePosix::file_write_meta",
// 	  "Trying to write metadata to the unopened file %s",
// 	  (path_ + name_).c_str());
//   }

//   // Determine the attribute rank

//   hsize_t meta_rank = 5;

//   if (n4 == 0) -- meta_rank;
//   if (n3 == 0) -- meta_rank;
//   if (n2 == 0) -- meta_rank;
//   if (n1 == 0) -- meta_rank;

//   ASSERT ("FilePosix::data_set","d is out of range",
// 	  (1 <= meta_rank && meta_rank <= 5));


//   hsize_t meta_size[5];

//   meta_size[0] = n0;
//   meta_size[1] = n1;
//   meta_size[2] = n2;
//   meta_size[3] = n3;
//   meta_size[4] = n4;
  
//   // Open the data space

//   hid_t meta_space_id = H5Screate_simple (meta_rank, meta_size, NULL);

//   hid_t meta_id = H5Acreate ( file_id_,
// 			      name.c_str(),
// 			      scalar_to_hdf5_(type),
// 			      meta_space_id,
// 			      H5P_DEFAULT);

//   H5Awrite (meta_id, scalar_to_hdf5_(type), buffer);

//   H5Sclose(meta_space_id);
//   H5Aclose(meta_id);
// }

// //----------------------------------------------------------------------

// void FilePosix::data_read_meta
//   ( void * buffer, std::string name,  enum scalar_type * type,
//     int * n0, int * n1, int * n2, int * n3, int * n4) throw()
// {
//   // error check file open

//   if ( ! is_file_open_) {

//     ERROR1("FilePosix::data_read_meta",
// 	  "Trying to read attribute from the unopened file %s",
// 	  (path_ + name_).c_str());
//   }

//   // error check dataset open

//   if (! is_data_open_) {

//     ERROR1("FilePosix::data_read_meta",
// 	  "Trying to read attribute from unopened dataset %s",
// 	  data_name_.c_str());
//   }

//   hid_t meta_id = H5Aopen_name(data_set_id_, name.c_str());

//   if (meta_id < 0) {

//     ERROR3("FilePosix::data_read_meta",
// 	  "H5Aopen_name() returned %d when opening attribute %s in file %s",
// 	  meta_id, name.c_str(),full_name.c_str());
//   }

//   // Get attribute size

//   hid_t meta_space_id = H5Aget_space (meta_id);

//   // error check rank

//   int rank = H5Sget_simple_extent_ndims(meta_space_id);

//   if (rank > 5) {

//     ERROR3("FilePosix::data_read_meta",
// 	  "Attribute %s in file %s has unsupported rank %d",
// 	  name.c_str(),name_.c_str(),rank);
//   }

//   // set output parameters

//   scalar_type scalar_type = hdf5_to_scalar_(H5Aget_type (meta_id));

//   if (type) (*type) = scalar_type;

//   hsize_t n[5];
//   H5Sget_simple_extent_dims(meta_space_id,n,0);

//   if (n0) (*n0) = n[0];
//   if (n1) (*n1) = n[1];
//   if (n2) (*n2) = n[2];
//   if (n3) (*n3) = n[3];
//   if (n4) (*n4) = n[4];

//   // Read the attribute

//   H5Aread(meta_id, scalar_to_hdf5_(scalar_type), buffer);

// }

// //----------------------------------------------------------------------

// void FilePosix::data_write_meta
//   ( const void * buffer, std::string name, enum scalar_type type,
//     int n0, int n1, int n2, int n3, int n4) throw()
// {
//   // error check file open

//   if ( ! is_file_open_) {

//     ERROR1("FilePosix::data_write_meta",
// 	  "Trying to write metadata to the unopened file %s",
// 	  (path_ + name_).c_str());
//   }

//   // error check dataset open

//   if (! is_data_open_) {

//     ERROR1("FilePosix::data_write_meta",
// 	  "Trying to read attribute from unopened dataset %s",
// 	  data_name_.c_str());
//   }

//   // Determine the attribute rank

//   hsize_t meta_rank = 5;

//   if (n4 == 0) -- meta_rank;
//   if (n3 == 0) -- meta_rank;
//   if (n2 == 0) -- meta_rank;
//   if (n1 == 0) -- meta_rank;

//   ASSERT ("FilePosix::data_set","d is out of range",
// 	  (1 <= meta_rank && meta_rank <= 5));


//   hsize_t meta_size[5];

//   meta_size[0] = n0;
//   meta_size[1] = n1;
//   meta_size[2] = n2;
//   meta_size[3] = n3;
//   meta_size[4] = n4;
  
//   // Open the data space

//   hid_t meta_space_id = H5Screate_simple (meta_rank, meta_size, NULL);

//   hid_t meta_id = H5Acreate ( data_set_id_,
// 			      name.c_str(),
// 			      scalar_to_hdf5_(type),
// 			      meta_space_id,
// 			      H5P_DEFAULT);

//   H5Awrite (meta_id, scalar_to_hdf5_(type), buffer);

//   H5Sclose(meta_space_id);
//   H5Aclose(meta_id);
// }  

// //----------------------------------------------------------------------

// void FilePosix::group_open (std::string name) throw()
// {
//   // Close current group if open
//   group_close();
  
//   // Error if not an absolute path name
//   if (name[0] != '/') {

//     ERROR1("FilePosix::data_write_meta",
// 	  "Group name '%s' must begin with '/'", name.c_str());
//   }
  
//   group_id_ = H5Gopen(file_id_, name.c_str());

//   is_group_open_ = true;
  
// }

// //----------------------------------------------------------------------

// void FilePosix::group_create (std::string group_path) throw()
// {
//   // Close current group if open
//   group_close();

//   // Error if not an absolute path name
//   if (group_path[0] != '/') {

//     ERROR1("FilePosix::data_write_meta",
// 	  "Group name '%s' must begin with '/'", group_path.c_str());
//   }


//   // Start at root group
  
//   std::string group_full = "/";
//   std::string group_rest = group_path;
//   group_rest.erase(0,1);

//   group_id_ = H5Gopen(file_id_,group_full.c_str());

//   // Loop through ancestor groups, creating if needed

//   size_t pos;

//   bool done = false;
//   while ( ! done ) {

//     pos = group_rest.find("/",0);

//     // Get the next subgroup name

//     std::string group = group_rest.substr(0,pos);
//     group_rest.erase(0,pos+1);

//     // Loop through children to find if subgroup exists

//     H5G_info_t group_info;
//     H5Gget_info_by_name (file_id_, group_full.c_str(), &group_info, H5P_DEFAULT);

//     bool group_exists = false;
//     for (int i=0; i<group_info.nlinks; i++) {
//       char group_name[80];
//       H5Lget_name_by_idx (file_id_, group_full.c_str(),
// 			  H5_INDEX_NAME,
// 			  H5_ITER_NATIVE,i,
// 			  group_name,80,H5P_DEFAULT);
//       if (group == group_name) {
// 	group_exists = true;
// 	break;
//       }
//     }

//     group_full = group_full + group + "/" ;

//     //  Open or create next group in path

//     hid_t group_new;

//     if (group_exists) {
//       group_new= H5Gopen(file_id_,group_full.c_str());
//     } else {
//       group_new = H5Gcreate (file_id_,group_full.c_str(), H5P_DEFAULT);
//     }

//     // Close parent group

//     H5Gclose (group_id_);

//     // Update group id

//     group_id_ = group_new;
    
//     done = (pos == std::string::npos);

//   }

//   is_group_open_ = true;
// }

// //----------------------------------------------------------------------

// void FilePosix::group_close () throw()
// {
//   herr_t status;
//   if (is_group_open_) {
//     //
//     status = H5Gclose(group_id_);
//   }
//   is_group_open_ = false;
// }

// //======================================================================

// int FilePosix::scalar_to_hdf5_ (enum scalar_type type) const throw()
// {
//   // (*) NATIVE    -   FLOAT DOUBLE LDOUBLE
//   // ( ) IEEE      -   F32BE F64BE     -
//   // ( ) STD     B16BE B32BE B64BE     -
//   // Types: http://www.hdfgroup.org/HDF5/Tutor/datatypes.html#native-types
//   // char          H5T_NATIVE_CHAR   H5T_STD_I8BE or H5T_STD_I8LE
//   // float         H5T_NATIVE_FLOAT   H5T_IEEE_F32BE or H5T_IEEE_F32LE  
//   // double        H5T_NATIVE_DOUBLE   H5T_IEEE_F64BE or H5T_IEEE_F64LE  
//   // unsigned char H5T_NATIVE_UCHAR   H5T_STD_U8BE or H5T_STD_U8LE
//   // int           H5T_NATIVE_INT   H5T_STD_I32BE or H5T_STD_I32LE
//   // short:        H5T_NATIVE_SHORT   H5T_STD_I16BE or H5T_STD_I16LE
//   // long:         H5T_NATIVE_LONG   H5T_STD_I32BE, H5T_STD_I32LE,
//   //               H5T_STD_I64BE or H5T_STD_I64LE
//   // long long:    H5T_NATIVE_LLONG   H5T_STD_I64BE or H5T_STD_I64LE

//   hid_t hdf5_type;
//   switch (type) {
//   case scalar_type_unknown:
//     ERROR("FilePosix::scalar_to_hdf5_",
// 	  "scalar_type_unknown not implemented");
//     hdf5_type = 0;
//     break;
//   case scalar_type_float:
//     hdf5_type = H5T_NATIVE_FLOAT;
//     break;
//   case scalar_type_double:
//     hdf5_type = H5T_NATIVE_DOUBLE;
//     break;
//   case scalar_type_long_double:
//     ERROR("FilePosix::scalartype_",
// p	  "long double not supported");
//     hdf5_type = 0;
//     break;
//   case scalar_type_char:
//     hdf5_type = H5T_NATIVE_CHAR;
//     //  H5T_STD_I8BE
//     //  H5T_STD_I8LE
//     break;
//   case scalar_type_int:
//     hdf5_type = H5T_NATIVE_INT;
//     //  H5T_STD_I32BE 
//     //  H5T_STD_I32LE
//     break;
//   case scalar_type_long:
//     //  H5T_STD_I64BE
//     //  H5T_STD_I64LE
//     hdf5_type = H5T_NATIVE_LONG;
//     break;
//   default:
//     ERROR1("FilePosix::type_", "unsupported type %d", type);
//     hdf5_type = 0;
//     break;
//   }
//   return hdf5_type;
// }

// //----------------------------------------------------------------------

// enum scalar_type FilePosix::hdf5_to_scalar_ (int hdf5_type) const throw()
// {

//   H5T_class_t hdf5_class = H5Tget_class(hdf5_type);
//   size_t      hdf5_size  = H5Tget_size (hdf5_type);

//   enum scalar_type type;

//   if (hdf5_class == H5T_INTEGER) {

//     if (hdf5_size == sizeof(char))   type = scalar_type_char;
//     if (hdf5_size == sizeof(int))    type = scalar_type_int;
//     if (hdf5_size == sizeof(long))   type = scalar_type_long;

//   } else if (hdf5_class == H5T_FLOAT) {

//     if (hdf5_size == sizeof(float))  type = scalar_type_float;
//     if (hdf5_size == sizeof(double)) type = scalar_type_double;

//   } else {

//     ERROR2("FilePosix::data_get",
// 	  "Unknown type of class %d and size %d",
// 	  hdf5_class, int(hdf5_size));
//   }

//   return type;

// }

// //----------------------------------------------------------------------

// size_t FilePosix::scalar_size_ (enum scalar_type type) const
// {
//   int bytes = 0;

//   switch (type) {
//   case scalar_type_unknown:
//     ERROR("FilePosix::scalar_size",
// 	  "scalar_type_unknown not implemented");
//     break;
//   case scalar_type_float:
//     bytes = sizeof (float);
//     break;
//   case scalar_type_double:
//     bytes = sizeof(double);
//     break;
//   case scalar_type_long_double:
//     bytes = sizeof(long double);
//     break;
//   case scalar_type_char:
//     bytes = sizeof(char);
//     break;
//   case scalar_type_int:
//     bytes = sizeof(int);
//     break;
//   case scalar_type_long:
//     bytes = sizeof(long);
//     break;
//   default:
//     ERROR1("FilePosix::scalar_size_",
// 	   "unsupported type %d", type);
//     hdf5_type = 0;
//     break;
//   }
//   return bytes;
// }

