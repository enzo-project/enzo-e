// $Id: disk_File.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_FILE_HPP
#define DISK_FILE_HPP

/// @file     disk_File.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the File class

class File {

  /// @class    File
  /// @ingroup  Disk
  /// @brief    Internal representation of a file on disk

public: // interface

  /// Constructor
  File() throw();

//   /// Destructor
//   ~File() throw();

//   /// Copy constructor
//   File(const File & diskfile) throw();

//   /// Assignment operator
//   File & operator= (const File & diskfile) throw();

  void open(std::string filename) throw();
  void close() throw();
  void read() throw();
  void write() throw();
  void flush() throw();

private: // functions


private: // attributes


};

#endif /* DISK_FILE_HPP */

