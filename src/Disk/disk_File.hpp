// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_FILE_HPP
#define DISK_FILE_HPP

/// @file     disk_File.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    [\ref Disk] Declaration of the File class

class File {

  /// @class    File
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Internal representation of a file on disk

public: // interface

  /// Constructor
  File() throw();

  /// Open the given named file
  void open(std::string filename) throw();

  /// Close the file
  void close() throw();

  /// Read data from the file
  void read() throw();

  /// Write data to the file
  void write() throw();

private: // functions


private: // attributes


};

#endif /* DISK_FILE_HPP */

