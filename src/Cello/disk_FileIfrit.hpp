// See LICENSE_CELLO file for license and copyright information

/// @file     disk_FileIfrit.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:05:34 PST 2008
/// @brief    [\ref Disk] Interface for the FileIfrit class
 
#ifndef DISK_IFRIT_HPP
#define DISK_IFRIT_HPP

class FileIfrit {

  /// @class    FileIfrit
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Class for writing and reading IFRIT files
  ///
  /// An FileIfrit object currently corresponds to a single IFRIT file /
  /// group dataset.  "IFrIT is a powerful tool that can be used to
  /// visualize 3-dimensional data sets."
  /// http://sites.google.com/site/ifrithome/

public: /// interface

  /// Initialize the FileIfrit object
  FileIfrit() {};

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

  }

  /// Read a 3d array from an ifrit file
  void read_bin  (std::string name, 
		  float *     buffer, 
		  int *       nx, 
		  int *       ny, 
		  int *       nz) throw ();

  /// Write a 3d array to an ifrit file
  void write_bin (std::string name, 
		  float *     buffer, 
		  int         nx, 
		  int         ny, 
		  int         nz) throw();

};

#endif /* DISK_IFRIT_HPP */

