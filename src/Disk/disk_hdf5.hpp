/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * SYNOPSIS:
 *
 *    
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */

#ifndef DISK_HDF5_HPP
#define DISK_HDF5_HPP

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
  void dataset_open_write (std::string name, 
			   int nx, int ny, int nz,
			   int ix0=0,int iy0=0,int iz0=0, 
			   int mx=0, int my=0, int mz=0);
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

