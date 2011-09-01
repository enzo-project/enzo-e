// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_FILE_HPP
#define DISK_FILE_HPP

/// @file     disk_File.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-26
/// @brief    [\ref Disk] Declaration of the File class

/// @enum scalar_type
/// @brief Simple scalar data type, e.g. scalar_int, scalar_float, etc.
enum scalar_type {
  scalar_type_unknown,
  scalar_type_char,        // Used for string data, with size + 1 for \0 
  scalar_type_int,
  scalar_type_long,
  scalar_type_float,
  scalar_type_double,
  scalar_type_long_double
};

#define MAX_DISK_ARRAY_RANK 5

class File {

  /// @class    File
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Internal representation of a file on disk

public: // interface

  /// Constor: create a file with the give path, file name, and access mode
  File (std::string path, std::string name, std::string mode) throw();

  /// Open the file
  virtual void file_open () throw() = 0;

  /// Close the file
  virtual void file_close () throw() = 0;
  
  // /// Set the current attribute type
  // virtual void attr_set
  // ( std::string name, 
  //   enum scalar_type type,
  //   int n0=1, int n1=1, int n2=1, int n3=1, int n4=1);
			      
  // /// Read attribute from the file
  // virtual void attr_read(void * buffer) throw() = 0;

  // /// Write attribute to the file
  // virtual void attr_write(const void * buffer) throw() = 0;

  /// Read data from the file
  virtual void data_read
  ( void * buffer, std::string name,  enum scalar_type * type,
    int * n0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw() = 0;

  /// Write data to the file
  virtual void data_write
  ( const void * buffer, std::string name, enum scalar_type type,
    int n0, int n1=1, int n2=1, int n3=1, int n4=1) throw() = 0;

protected: // attributes

  /// Path to the file
  std::string path_;

  /// Name of the file
  std::string name_;

  /// Intended mode for the file, "r", "w", or "rw"
  std::string mode_;

};

#endif /* DISK_FILE_HPP */

