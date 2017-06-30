// See LICENSE_CELLO file for license and copyright information

/// @file     disk_File.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-26
/// @brief    [\ref Disk] Declaration of the File class

#ifndef DISK_FILE_HPP
#define DISK_FILE_HPP

#define MAX_DISK_ARRAY_RANK 5

class File {

  /// @class    File
  /// @ingroup  Disk
  /// @brief    [\ref Disk] Internal representation of a file on disk

public: // interface

  /// Create a file with the given path and filename
  File (std::string path, std::string name) throw()
    : path_(path),
      name_(name)
  {}

  /// Create a file with the given path and filename
  virtual ~File () throw()
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | path_;
    p | name_;
  }

public: // virtual functions

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
  ( void * buffer, std::string name,  int * s_type,
    int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw() = 0;
  
  /// Write a metadata item associated with the file
  virtual void file_write_meta
  ( const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw() = 0;
  

  // Datasets

  /// Create a new dataset for writing (and open it)
  virtual void data_create
  ( std::string name,  int type,
    int m1=1, int m2=0, int m3=0, int m4=0,
    int n1=0, int n2=0, int n3=0, int n4=0,
    int o1=0, int o2=0, int o3=0, int o4=0) throw() = 0;

  /// Open an existing dataset for reading
  virtual void data_open
  ( std::string name,  int * type,
    int * m1=0, int * m2=0, int * m3=0, int * m4=0) throw() = 0;

  /// Select a subset of the data
  virtual void data_slice
  ( int m1, int m2, int m3, int m4,
    int n1, int n2, int n3, int n4,
    int o1, int o2, int o3, int o4) throw() = 0;

  /// Return the size of the disk dataset
  virtual int data_size (int * m4) throw() = 0;

  /// Read from the opened dataset
  virtual void data_read (void * buffer) throw() = 0;

  /// Write to the opened dataset
  virtual void data_write 
  (const void * buffer) throw() = 0;

  /// Close the opened dataset
  virtual void data_close () throw() = 0;

  /// Read a metadata item associated with the opened dataset
  virtual void data_read_meta
  ( void * buffer, std::string name,  int * s_type,
    int * n1, int * n2=0, int * n3=0, int * n4=0) throw() = 0;
  
  /// Write a metadata item associated with the opened dataset
  virtual void data_write_meta
  ( const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw() = 0;

  /// Create memory space
  virtual void mem_create
  ( int mx, int my, int mz,
    int nx, int ny, int nz,
    int gx, int gy, int gz ) = 0;

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
  ( void * buffer, std::string name,  int * s_type,
    int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw() = 0;
  
  /// Write a metadata item associated with the opened group
  virtual void group_write_meta
  ( const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw() = 0;

protected: // attributes

  /// Path to the file
  std::string path_;

  /// Name of the file
  std::string name_;

};

#endif /* DISK_FILE_HPP */

