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

  /// Create an HDF5 file with the given path and filename

  FileHdf5 (std::string path, std::string name) throw();

  //--------------------------------------------------
  // Files
  //--------------------------------------------------

  /// Open an existing file

  virtual void file_open () throw();

  /// Create a new file

  virtual void file_create () throw();

  /// Close the file

  virtual void file_close () throw();
  
  /// Read a metadata item associated with the file

  virtual void file_meta_read
  ( void * buffer, std::string name,  enum scalar_type * s_type,
    int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();
  
  /// Write a metadata item associated with the file

  virtual void file_meta_write
  ( const void * buffer, std::string name, enum scalar_type type,
    int n0=1, int n1=0, int n2=0, int n3=0, int n4=0) throw();
  
  //--------------------------------------------------
  // Datasets
  //--------------------------------------------------

  /// Open an existing dataset for reading

  virtual void data_open
  ( std::string name,  enum scalar_type * type,
    int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  /// Create a new dataset for writing (and open it)

  virtual void data_create
  ( std::string name,  enum scalar_type type,
    int n0=1, int n1=0, int n2=0, int n3=0, int n4=0) throw();

  /// Read from the opened dataset

  virtual void data_read (void * buffer) throw();

  /// Write to the opened dataset

  virtual void data_write (const void * buffer) throw();

  /// Close the opened dataset

  virtual void data_close () throw();

  /// Read a metadata item associated with the opened dataset

  virtual void data_meta_read
  ( void * buffer, std::string name,  enum scalar_type * s_type,
    int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();
  
  /// Write a metadata item associated with the opened dataset

  virtual void data_meta_write
  ( const void * buffer, std::string name, enum scalar_type type,
    int n0=1, int n1=0, int n2=0, int n3=0, int n4=0) throw();
  
private: // functions

  /// Convert the scalar type to HDF5 datatype
  int scalar_to_hdf5_(enum scalar_type type) const throw();

  /// Convert the scalar type to an HDF5 datatype
  enum scalar_type hdf5_to_scalar_(int type) const throw();

private: // attributes

  /// HDF5 file descriptor
  hid_t file_id_;

  /// HDF5 dataset descriptor
  hid_t data_set_id_;

  /// HDF5 dataspace descriptor
  hid_t data_space_id_;

  /// HDF5 attribute descriptor
  hid_t attribute_id_;

  /// HDF5 error satus
  herr_t status_id_;

  /// HDF5 dataset name
  std::string data_name_;

  /// Type of data in the HDF5 datatype
  scalar_type data_type_;

  /// Dataset rank, 0 to 5
  int data_rank_;

  /// Dataset size
  hsize_t data_size_[5];

  /// Whether file is open or closed
  bool  is_file_open_;

  /// Whether dataset is open or closed
  bool  is_data_open_;

};

#endif /* DISK_FILE_HDF5_HPP */

