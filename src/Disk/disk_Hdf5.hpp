// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_HDF5_HPP
#define DISK_HDF5_HPP

/// @file     disk_hdf5.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:05:34 PST 2008
/// @todo     Refactor interface to be hdf5-independent (groups, datasets, etc.)
/// @todo     Rename to, e.g., FileHdf5
/// @todo     Add support for relative/absolute directory / group
/// @todo     Add unit tests for group operations to test_disk_hdf5
/// @brief    Interface for the Hdf5 class

class Hdf5 {

  /// @class    Hdf5
  /// @ingroup  Disk
  /// @brief    Class for writing and reading HDF5 files
  ///
  /// An Hdf5 object currently corresponds to a single HDF5 file / group

public: // interface

  /// Initialize the Hdf5 object
  Hdf5();

  /// Open the file with the given mode
  int file_open  (std::string name, std::string mode);

  /// Close the file
  void file_close ();

  /// Open the given group
  void group_open (std::string name);

  /// Close the current group
  void group_close ();

  /// Open the given dataset with given size for reading
  void dataset_open_read (std::string name, int * nx, int * ny, int * nz);

  /// Open the given dataset with the given size for writing
  void dataset_open_write (std::string name, 
			   int nx, int ny, int nz);

  /// Close the current dataset
  void dataset_close ();

  /// Read the current dataset into the buffer
  void read  (Scalar * buffer);

  /// Write the current dataset from the buffer
  void write (Scalar * buffer);

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

  /// HDF5 data type
  hid_t datatype_;

  /// Last error
  herr_t      status_;

};

#endif /* DISK_HDF5_HPP */

