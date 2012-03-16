// See LICENSE_CELLO file for license and copyright information

/// @file     disk_FileHdf5.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:05:34 PST 2008
/// @brief    [\ref Disk] Interface for the FileHdf5 class

#ifndef DISK_FILE_HDF5_HPP
#define DISK_FILE_HDF5_HPP

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
  ~FileHdf5 () throw();

  // Files

  /// Open an existing file
  virtual void file_open () throw();

  /// Create a new file
  virtual void file_create () throw();

  /// Close the file
  virtual void file_close () throw();
  
  /// Read a metadata item associated with the file
  virtual void file_read_meta
  ( void * buffer, std::string name,  enum scalar_type * s_type,
    int * nx=0, int * ny=0, int * nz=0) throw();
  
  /// Write a metadata item associated with the file

  void file_write_meta
  ( const void * buffer, std::string name, enum scalar_type type,
    int nx=1, int ny=0, int nz=0) throw()
  {  write_meta_ ( file_id_, buffer, name, type, nx,ny,nz); }

  // Datasets

  /// Open an existing dataset for reading
  virtual void data_open
  ( std::string name,  enum scalar_type * type,
    int * nx=0, int * ny=0, int * nz=0) throw();

  /// Create a new dataset for writing (and open it)
  virtual void data_create
  ( std::string name,  enum scalar_type type,
    int nx=1, int ny=0, int nz=0) throw();

  /// Read from the opened dataset
  virtual void data_read (void * buffer) throw();

  /// Write to the opened dataset
  virtual void data_write (const void * buffer) throw();

  /// Close the opened dataset
  virtual void data_close () throw();

  /// Read a metadata item associated with the opened dataset
  virtual void data_read_meta
  ( void * buffer, std::string name,  enum scalar_type * s_type,
    int * nx=0, int * ny=0, int * nz=0) throw();
  
  /// Write a metadata item associated with the opened dataset
  virtual void data_write_meta
  ( const void * buffer, std::string name, enum scalar_type type,
    int nx=1, int ny=0, int nz=0) throw();


  // Groups

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
  ( void * buffer, std::string name,  enum scalar_type * s_type,
    int * nx=0, int * ny=0, int * nz=0) throw();
  
  /// Write a metadata item associated with the opened group
  void group_write_meta
  ( const void * buffer, std::string name, enum scalar_type type,
    int nx=1, int ny=0, int nz=0) throw()
  {  write_meta_ ( group_id_, buffer, name, type, nx,ny,nz ); }

  /// Set the compression level
  void set_compress (int level) throw ();

  /// Return the compression level
  int compress () throw () {return compress_level_; }


protected: // functions

  virtual void write_meta_
  ( hid_t id, const void * buffer, std::string name, enum scalar_type type,
    int nx=1, int ny=0, int nz=0) throw();

private: // functions

  /// Convert the scalar type to HDF5 datatype
  int scalar_to_hdf5_(enum scalar_type type) const throw();

  /// Convert the scalar type to an HDF5 datatype
  enum scalar_type hdf5_to_scalar_(int type) const throw();

  /// Convert a relative path to an absolute path
  std::string relative_to_absolute_
  (std::string path_relative, std::string path_current) const throw();

  /// Get output extents
  void get_output_extents_
  ( hid_t space_id, int * nx, int * ny, int * nz) throw();

  /// Determine rank
  hid_t create_dataspace_(int nx,int ny,int nz, hsize_t * data_size) throw();

  /// Close the given dataspace
  void close_dataspace_ (hid_t space_id) throw();

  /// Open the dataset
  hid_t open_dataset_ (hid_t group, std::string name) throw();

  /// Close the dataset
  void close_dataset_ () throw();

  /// Return the space for the given dataset
  hid_t get_data_space_(hid_t dataset_id, std::string name) throw ();

  /// Return the space for the given attribute
  hid_t get_attr_space_(hid_t dataset_id, std::string name) throw ();

private: // attributes

  /// HDF5 file descriptor
  hid_t file_id_;

  /// Whether file is open or closed
  bool  is_file_open_;


  /// HDF5 dataset descriptor
  hid_t data_id_;

  /// HDF5 dataspace descriptor
  hid_t space_id_;


  /// HDF5 attribute descriptor
  hid_t attribute_id_;


  /// HDF5 group descriptor
  hid_t group_id_;

  /// Group name 
  std::string group_name_;

  /// HDF5 group property list
  hid_t group_prop_;

  /// Whether a group is open or closed
  bool is_group_open_;

  /// HDF5 dataset name
  std::string data_name_;

  /// Type of data in the HDF5 datatype
  scalar_type data_type_;

  /// Dataset rank, 0 to 5
  int data_rank_;

  /// Dataset size
  hsize_t data_size_[5];

  /// HDF5 dataset property list
  hid_t data_prop_;

  /// Whether a dataset is open or closed
  bool  is_data_open_;

  /// Compression level
  int compress_level_;

};

#endif /* DISK_FILE_HDF5_HPP */

