#ifndef DISK_HDF5_HPP
#define DISK_HDF5_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      hdf5.hpp
 * @brief     Definition of the Hdf5 class
 * @author    James Bordner
 * @date      Thu Feb 21 16:05:34 PST 2008
 *
 * $Id$
 *
 *********************************************************************
 */
 
class Hdf5 {

/** 
 *********************************************************************
 *
 * @class     Hdf5
 * @brief     Class for writing and reading HDF5 files
 * @ingroup   Storage
 *
 * An Hdf5 object currently corresponds to a single HDF5 file / group
 * / dataset.
 *
 *********************************************************************
 */

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

public:

  /// Initialize the Hdf5 object
  Hdf5();
  int file_open  (std::string name, std::string mode);
  void file_close ();
  void group_open (std::string name);
  void group_close ();
  void dataset_open_read (std::string name, int * nx, int * ny, int * nz);
  void dataset_open_write (std::string name, int nx, int ny, int nz);
  void dataset_set_hyperslab (int ix0,int iy0, int iz0, 
			       int mx, int my, int mz);
  void dataset_close ();
  void read  (Scalar * buffer);
  void write (Scalar * buffer);

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

private:


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

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

private:

};

#endif /* DISK_HDF5_HPP */

