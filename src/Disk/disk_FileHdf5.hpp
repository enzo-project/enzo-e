// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_FILE_HDF5_HPP
#define DISK_FILE_HDF5_HPP

/// @file     disk_FileHdf5.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:05:34 PST 2008
/// @todo     Refactor interface to be hdf5-independent (groups, datasets, etc.)
/// @todo     Add support for [relative|absolute] [directories|groups] (?)
/// @todo     Support multiple float types: std, native, ieee
/// @todo     Add error handling (see H5E API)
/// @todo     Add state checks for file open before adding dataset, etc.
/// @todo     Add support for compression (see H5Z API)
/// @brief    [\ref Disk] Interface for the FileHdf5 class

class FileHdf5 {

  /// @class    FileHdf5
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Class for writing and reading HDF5 files
  ///
  /// An FileHdf5 object currently corresponds to a single HDF5 file / group

public: // interface

  /// Initialize the FileHdf5 object
  FileHdf5();

  /// Open the file with the given mode
  int open_file  (std::string name, std::string mode);

  /// Close the file
  void close_file ();

  /// Open the given group
  void open_group (std::string name);

  /// Close the current group
  void close_group ();

  /// Open the given dataset with given size for reading
  void open_dataset (std::string name, 
		     int * nx, int * ny, int * nz);

  /// Open the given dataset with the given size for writing
  void open_dataset (std::string name, 
		     enum precision_enum precision,
		     int nx, int ny, int nz);

  /// Close the current dataset
  void close_dataset ();

  /// Read the current dataset into the buffer
  void read  (char              * buffer,
	      enum precision_enum precision);

  /// Write the current dataset from the buffer
  void write (char              * buffer,
	      enum precision_enum precision);

private: // functions

  /// Return the HDF5 datatype for the given precision
  int datatype_(enum precision_enum precision);

private: // attributes

  /// HDF5 file descriptor
  hid_t file_;

  /// HDF5 file name
  std::string file_name_;

  /// HDF5 file mode
  std::string file_mode_;

  /// Whether file is open or closed
  bool  is_file_open_;

  /// HDF5 dataset descriptor
  hid_t dataset_;

  /// HDF5 dataset name
  std::string dataset_name_;

  /// HDF5 dataspace descriptor
  hid_t dataspace_;

  /// Last error
  herr_t      status_;

};

#endif /* DISK_FILE_HDF5_HPP */

