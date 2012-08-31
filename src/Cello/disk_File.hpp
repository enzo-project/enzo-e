// See LICENSE_CELLO file for license and copyright information

/// @file     disk_File.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-26
/// @brief    [\ref Disk] Declaration of the File class

#ifndef DISK_FILE_HPP
#define DISK_FILE_HPP

/// @enum scalar_enum
/// @brief Simple scalar data type, e.g. scalar_int, scalar_float, etc.
enum scalar_enum {
  scalar_type_unknown,
  scalar_type_char,        // Used for string data, with size + 1 for \0 
  scalar_type_int,
  scalar_type_long,
  scalar_type_long_long,
  scalar_type_float,
  scalar_type_double,
  scalar_type_long_double
};

typedef char scalar_type;

#define MAX_DISK_ARRAY_RANK 5

class File {

  /// @class    File
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Internal representation of a file on disk

public: // interface

  /// Create a file with the given path and filename
  File (std::string path, std::string name) throw();

  /// Create a file with the given path and filename
  virtual ~File () throw()
  {}

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    p | path_;
    p | name_;
  }
#endif

  //--------------------------------------------------
  // Files
  //--------------------------------------------------

  /// Open an existing file
  virtual void file_open () throw() = 0;

  /// Create a new file
  virtual void file_create () throw() = 0;

  /// Close the file
  virtual void file_close () throw() = 0;
  
  /// Read a metadata item associated with the file
  virtual void file_read_meta
  ( void * buffer, std::string name,  scalar_type * s_type,
    int * nx=0, int * ny=0, int * nz=0) throw() = 0;
  
  /// Write a metadata item associated with the file
  virtual void file_write_meta
  ( const void * buffer, std::string name, scalar_type type,
    int nx=1, int ny=0, int nz=0) throw() = 0;
  

  // Datasets

  /// Open an existing dataset for reading
  virtual void data_open
  ( std::string name,  scalar_type * type,
    int * nx=0, int * ny=0, int * nz=0) throw() = 0;

  /// Create a new dataset for writing (and open it)
  virtual void data_create
  ( std::string name,  scalar_type type,
    int nxd=1, int nyd=0, int nzd=0,
    int nx=0,  int ny=0,  int nz=0) throw() = 0;

  /// Read from the opened dataset
  virtual void data_read (void * buffer) throw() = 0;

  /// Write to the opened dataset
  virtual void data_write 
  (const void * buffer) throw() = 0;

  /// Close the opened dataset
  virtual void data_close () throw() = 0;

  /// Read a metadata item associated with the opened dataset
  virtual void data_read_meta
  ( void * buffer, std::string name,  scalar_type * s_type,
    int * nx, int * ny=0, int * nz=0) throw() = 0;
  
  /// Write a metadata item associated with the opened dataset
  virtual void data_write_meta
  ( const void * buffer, std::string name, scalar_type type,
    int nx=1, int ny=0, int nz=0) throw() = 0;


  // Groups

  /// Return the number of subgroups in the current group
  virtual int group_count () const throw() = 0;

  /// Return the name of the ith subgroup
  virtual std::string group_name (size_t i) const throw() = 0;

  /// Change group name for subsequent open or create
  virtual void group_chdir (std::string name) throw() = 0;

  /// Open the existing group from current group name
  virtual void group_open () throw() = 0;

  /// Create a new group from current group name (and open it)
  virtual void group_create () throw() = 0;

  /// Get the current group
  virtual void group_close () throw() = 0;

  /// Read a metadata item associated with the opened group
  virtual void group_read_meta
  ( void * buffer, std::string name,  scalar_type * s_type,
    int * nx, int * ny=0, int * nz=0) throw() = 0;
  
  /// Write a metadata item associated with the opened group
  virtual void group_write_meta
  ( const void * buffer, std::string name, scalar_type type,
    int nx=1, int ny=0, int nz=0) throw() = 0;

protected: // attributes

  /// Path to the file
  std::string path_;

  /// Name of the file
  std::string name_;

};

#endif /* DISK_FILE_HPP */

