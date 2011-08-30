// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_FILE_HPP
#define DISK_FILE_HPP

/// @file     disk_File.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    [\ref Disk] Declaration of the File class

/// @enum scalar_type
/// @brief Simple scalar data type, e.g. scalar_int, scalar_float, etc.
enum scalar_type {
  scalar_char,        // Used for string data, with size + 1 for \0 
  scalar_int,
  scalar_long_int,
  scalar_float,
  scalar_double,
  scalar_long_double
};

class File {

  /// @class    File
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Internal representation of a file on disk

public: // interface

  /// Constor: create a file with the give path, file name, and access mode
  File(std::string path, std::string name, std::string mode) throw();

  /// Open the file
  virtual int open() throw() = 0;

  /// Close the file
  virtual void close() throw() = 0;
  
  /// Set the current attribute type
  virtual void set_attr_type 
  ( enum scalar_type scalar,
    int n0=1, int n1=1, int n2=1, int n3=1, int n4=1);
			      
  virtual void set_data_type 
  ( enum scalar_type scalar,
    int n0=1, int n1=1, int n2=1, int n3=1, int n4=1);
			      
  /// Read data from the file
  virtual void read_data(char * buffer) throw() = 0;

  /// Write data to the file
  virtual void write_data(const char * buffer) throw() = 0;

  /// Read attribute from the file
  virtual void read_attr(char * buffer) throw() = 0;

  /// Write attribute to the file
  virtual void write_attr(const char * buffer) throw() = 0;

protected: // functions

  std::string path_;
  std::string name_;
  std::string mode_;

private: // attributes


};

#endif /* DISK_FILE_HPP */

