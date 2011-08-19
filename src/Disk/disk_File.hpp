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
  virtual int open_file(std::string filename, std::string mode) throw() = 0;

  /// Close the file
  virtual void close_file() throw() = 0;
  
  /// Read data from the file
  virtual void read(char * buffer, enum precision_enum precision) throw() = 0;

  /// Write data to the file
  virtual void write(const char * buffer, enum precision_enum precision) throw() = 0;

private: // functions


private: // attributes


};

#endif /* DISK_FILE_HPP */

