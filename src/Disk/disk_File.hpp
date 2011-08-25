// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_FILE_HPP
#define DISK_FILE_HPP

/// @file     disk_File.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    [\ref Disk] Declaration of the File class

/// @enum  file_content_type
/// @brief Argument for read() and write() (in Hierarchy, Patch, Block,
/// etc.) to specify what to read or write.

enum file_content_type {
  file_content_header,
  file_content_data,
  file_content_all
};

class File {

  /// @class    File
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Internal representation of a file on disk

public: // interface

  /// Constructor
  File() throw();

  /// Open the given named file
  virtual int open(std::string filename, std::string mode) throw() = 0;

  /// Close the file
  virtual void close() throw() = 0;
  
  /// Read data from the file
  virtual void read(char * buffer, enum precision_enum precision) throw() = 0;

  /// Write data to the file
  virtual void write(const char * buffer, enum precision_enum precision) throw() = 0;

private: // functions


private: // attributes


};

#endif /* DISK_FILE_HPP */

