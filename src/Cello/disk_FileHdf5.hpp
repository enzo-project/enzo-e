// See LICENSE_CELLO file for license and copyright information

/// @file     disk_FileHdf5.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:05:34 PST 2008
/// @brief    [\ref Disk] Interface for the FileHdf5 class

#ifndef DISK_FILE_HDF5_HPP
#define DISK_FILE_HDF5_HPP

//----------------------------------------------------------------------
/// @typedef hdf5_id
/// @brief  Define a type to hold an HDF5 "hid_t" to localize dependencies

typedef int64_t hdf5_id;

class FileHdf5 : public File {

  /// @class    FileHdf5
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Class for writing and reading HDF5 files
  ///
  /// An FileHdf5 object currently corresponds to a single HDF5 file / group

public: // interface

  /// Create an HDF5 file with the given path and filename

  FileHdf5 (std::string path, std::string name) throw();

  /// Destructor
  virtual ~FileHdf5 () throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    File::pup(p);

    p | file_id_;
    p | is_file_open_;
    p | data_id_;
    p | data_space_id_;
    p | mem_space_id_;
    p | attribute_id_;
    p | group_id_;
    p | group_name_;
    p | group_prop_;
    p | is_group_open_;
    p | data_name_;
    p | data_type_;
    p | data_rank_;
    p | data_prop_;
    p | is_data_open_;
    p | compress_level_;
  }

public: // virtual functions

  // Files

  /// Open an existing file
  virtual void file_open () throw();

  /// Create a new file
  virtual void file_create () throw();

  /// Close the file
  virtual void file_close () throw();
  
  /// Read scalar metadata item associated with the file
  virtual void file_read_scalar
  ( void * buffer, std::string name,  int * s_type) throw();

  /// Read a metadata item associated with the file
  virtual void file_read_meta
  ( void * buffer, std::string name,  int * s_type,
    int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();
  
  /// Write a metadata item associated with the file
  virtual void file_write_meta
  ( const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw()
  { write_meta_ ( file_id_, buffer, name, type, n1,n2,n3,n4); }

  // Datasets

/// Create a new dataset for writing (and open it)
  virtual void data_create
  ( std::string name,  int type,
    int m1=1, int m2=0, int m3=0, int m4=0,
    int n1=0, int n2=0, int n3=0, int n4=0,
    int o1=0, int o2=0, int o3=0, int o4=0) throw();

  /// Open an existing dataset for reading
  virtual void data_open
  ( std::string name,  int * type,
    int * m1=0, int * m2=0, int * m3=0, int * m4=0) throw();

  /// Select a subset of the data
  virtual void data_slice
  ( int m1, int m2, int m3, int m4,
    int n1, int n2, int n3, int n4,
    int o1, int o2, int o3, int o4) throw();

  /// Return the size of the disk dataset
  virtual int data_size (int * m4) throw();
  
  /// Read from the opened dataset
  virtual void data_read (void * buffer) throw();

  /// Write to the opened dataset
  virtual void data_write 
  (const void * buffer) throw();

  /// Close the opened dataset
  virtual void data_close () throw();

  /// Read a metadata item associated with the opened dataset
  virtual void data_read_meta
  ( void * buffer, std::string name,  int * s_type,
    int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();
  
  /// Write a metadata item associated with the opened dataset
  virtual void data_write_meta
  ( const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw()
  { write_meta_ ( data_id_, buffer, name, type, n1,n2,n3,n4); }

  /// Create memory space
  virtual void mem_create
  ( int mx, int my, int mz,
    int nx, int ny, int nz,
    int gx, int gy, int gz );

  virtual void mem_close ()
  { space_close_(mem_space_id_) ; }
  
  // Groups

  /// Number of subgroups in the current group
  virtual int group_count () const throw();

  /// Return the name of the ith subgroup
  virtual std::string group_name (size_t i) const throw();

  /// Change to the named group
  virtual void group_chdir (std::string name) throw();

  /// Open an existing group named buy group_chdir()
  virtual void group_open () throw();

  /// Create a new group named by group_chdir() (and open it)
  virtual void group_create () throw();

  /// Get the current group
  virtual void group_close () throw();
  
  /// Read a metadata item associated with the opened group
  void group_read_meta
  ( void * buffer, std::string name,  int * s_type,
    int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();
  
  /// Write a metadata item associated with the opened group
  virtual void group_write_meta
  ( const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw()
  { write_meta_ ( group_id_, buffer, name, type, n1,n2,n3,n4); }

public: // functions

  /// Set the compression level
  void set_compress (int level) throw ();

  /// Return the compression level
  int compress () throw () {return compress_level_; }

  /// Allocate a buffer for reading in a dataset of the given
  /// length and type
  char * allocate_buffer (int n, int type_data)
  {
    char * data;
    if (type_data == type_single) {
      data = (char *)new float [n];
    } else if (type_data == type_double) {
      data = (char *)new double [n];
    } else if (type_data == type_int8) {
      data = (char *)new int8_t [n];
    } else if (type_data == type_int16) {
      data = (char *)new int16_t [n];
    } else if (type_data == type_int32) {
      data = (char *)new int32_t [n];
    } else if (type_data == type_int64) {
      data = (char *)new int64_t [n];
    } else {
      data = nullptr;
      ERROR1 ("EnzoInitialHdf5::allocate_array_()",
              "Unsupported data type %d",type_data);
    }
    return data;
  }
  
protected: // functions

  virtual void write_meta_
  ( hdf5_id id, const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw();

private: // functions

  /// Convert the scalar type to HDF5 datatype
  hdf5_id scalar_to_hdf5_(int type) const throw();

  /// Convert the scalar type to an HDF5 datatype
  int hdf5_to_scalar_(hdf5_id type) const throw();

  /// Convert a relative path to an absolute path
  std::string relative_to_absolute_
  (std::string path_relative, std::string path_current) const throw();

  /// Get output extents
  void get_extents_
  ( hdf5_id space_id, int * n1, int * n2=0, int * n3=0, int * n4=0) throw();

  /// create the space for the array on disk

  /// create data spaces for memory or disk data
  hdf5_id space_create_
  (int m1, int m2, int m3, int m4,
   int n1, int n2, int n3, int n4,
   int o1, int o2, int o3, int o4) throw();

  /// given a space, select a slice
  hdf5_id space_slice_
  (hdf5_id space_id,
   int m1, int m2, int m3, int m4,
   int n1, int n2, int n3, int n4,
   int o1, int o2, int o3, int o4) throw();
  
  /// Close the given dataspace
  void space_close_ (hdf5_id space_id) throw();

  /// Return the space for the given dataset
  hdf5_id get_data_space_(hdf5_id dataset_id, std::string name) throw ();

  /// Close the given memspace
  void close_mem_space_ (hdf5_id space_id) throw();

  /// Return the space for the given attribute
  hdf5_id get_attr_space_(hdf5_id dataset_id, std::string name) throw ();

  /// Open the dataset
  hdf5_id open_dataset_ (hdf5_id group, std::string name) throw();

  /// Close the dataset
  void close_dataset_ () throw();

public: // static attributes

  /// Nodal list of files opened
  static std::map<const std::string,FileHdf5 *> file_list;

private: // attributes

  /// HDF5 file descriptor
  hdf5_id file_id_;

  /// Whether file is open or closed
  bool  is_file_open_;


  /// HDF5 dataset descriptor
  hdf5_id data_id_;

  /// HDF5 dataspace descriptor
  hdf5_id data_space_id_;

  /// HDF5 memory space descriptor
  hdf5_id mem_space_id_;


  /// HDF5 attribute descriptor
  hdf5_id attribute_id_;


  /// HDF5 group descriptor
  hdf5_id group_id_;

  /// Group name 
  std::string group_name_;

  /// HDF5 group property list
  hdf5_id group_prop_;

  /// Whether a group is open or closed
  bool is_group_open_;

  /// HDF5 dataset name
  std::string data_name_;

  /// Type of data in the HDF5 datatype
  int data_type_;

  /// Dataset rank, 1 to 4
  int data_rank_;

  /// HDF5 dataset property list
  hdf5_id data_prop_;

  /// Whether a dataset is open or closed
  bool  is_data_open_;

  /// Compression level
  int compress_level_;

};

#endif /* DISK_FILE_HDF5_HPP */

