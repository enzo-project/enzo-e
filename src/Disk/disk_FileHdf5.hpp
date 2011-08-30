// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_FILE_HDF5_HPP
#define DISK_FILE_HDF5_HPP

/// @file     disk_FileHdf5.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:05:34 PST 2008
/// @todo     Support multiple float types: std, native, ieee
/// @todo     Add error handling (see H5E API)
/// @todo     Add support for compression (see H5Z API)
/// @brief    [\ref Disk] Interface for the FileHdf5 class

class FileHdf5 : public File {

  /// @class    FileHdf5
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Class for writing and reading HDF5 files
  ///
  /// An FileHdf5 object currently corresponds to a single HDF5 file / group

public: // interface

  /// Initialize the FileHdf5 object
  FileHdf5(std::string path, std::string name, std::string mode) throw();

  /// Open the file with the given mode
  virtual void open () throw();

  /// Close the file
  virtual void close () throw();

  // /// Set the current attribute type
  // virtual void set_attr_type 
  // ( enum scalar_type scalar,
  //   int n0=1, int n1=1, int n2=1, int n3=1, int n4=1);
			      
  /// Read attribute from the file
  // virtual void read_attr(void * buffer) throw() = 0;

  /// Write attribute to the file
  //  virtual void write_attr(const void * buffer) throw() = 0;

  virtual void data_get
  ( std::string name, enum scalar_type * type, 
    int * n0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  virtual void data_set 
  ( std::string name, enum scalar_type type, 
    int n0, int n1=0, int n2=0, int n3=0, int n4=0) throw();

  /// Read data from the file
  virtual void data_read (void * buffer) throw() = 0;

  /// Write data to the file
  virtual void data_write (const void * buffer) throw() = 0;

  // /// Open the given group
  // void open_group (std::string group) throw();

  // /// Close the current group
  // void close_group () throw();

private: // functions

  /// Convert the scalar type to HDF5 datatype
  int type_(enum scalar_type type) throw();

private: // attributes

  /// HDF5 file descriptor
  hid_t file_id_;

  /// HDF5 dataset descriptor
  hid_t data_set_id_;

  /// HDF5 dataspace descriptor
  hid_t data_space_id_;

  /// HDF5 error satus
  herr_t status_id_;

  /// Whether file is open or closed
  bool  is_file_open_;

  /// HDF5 dataset name
  std::string data_name_;

  /// Type of data in the HDF5 datatype
  scalar_type data_type_;

  /// Dataset rank, 0 to 5
  int data_rank_;

  /// Dataset size
  hsize_t data_size_[5];

};

#endif /* DISK_FILE_HDF5_HPP */

