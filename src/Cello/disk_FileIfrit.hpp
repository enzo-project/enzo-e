// See LICENSE_CELLO file for license and copyright information

/// @file     disk_FileIfrit.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:05:34 PST 2008
/// @brief    [\ref Disk] Interface for the FileIfrit class
 
#ifndef DISK_IFRIT_HPP
#define DISK_IFRIT_HPP

class FileIfrit : public File {

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
  FileIfrit(std::string path, std::string name) throw();

  /// Destructor
  virtual ~FileIfrit () throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

  }

public: // virtual functions

  //--------------------------------------------------
  // Files
  //--------------------------------------------------

  /// Open an existing file
  virtual void file_open () throw();

  /// Create a new file
  virtual void file_create () throw();

  /// Close the file
  virtual void file_close () throw();
  
  /// Read a metadata item associated with the file
  virtual void file_read_meta
  ( void * buffer, std::string name,  int * s_type,
    int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();
  
  /// Write a metadata item associated with the file
  virtual void file_write_meta
  ( const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw();
  

  // Datasets

  /// Create a new dataset for writing (and open it)
  virtual void data_create
  ( std::string name,  int type,
    int m1=1, int m2=0, int m3=0, int m4=0,
    int n1=0, int n2=0, int n3=0, int n4=0,
    int o1=0, int o2=0, int o3=0, int o4=0) throw();

  /// Open an existing dataset for reading
  virtual void data_open
  ( std::string name,  int * type,
    int * m1=0, int * m2=0, int * m3=0, int * m4=0) throw();

  /// Select a subset of the data
  virtual void data_slice
  ( int m1, int m2, int m3, int m4,
    int n1, int n2, int n3, int n4,
    int o1, int o2, int o3, int o4) throw();

  /// Read from the opened dataset
  virtual void data_read (void * buffer) throw();

  /// Return the size of the disk dataset
  int data_size (int * m4) throw()
  { return 0; }
  
  /// Write to the opened dataset
  virtual void data_write 
  (const void * buffer) throw();

  /// Close the opened dataset
  virtual void data_close () throw();

  /// Read a metadata item associated with the opened dataset
  virtual void data_read_meta
  ( void * buffer, std::string name,  int * s_type,
    int * n1, int * n2=0, int * n3=0, int * n4=0) throw();
  
  /// Write a metadata item associated with the opened dataset
  virtual void data_write_meta
  ( const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw();

  /// Create memory space
  virtual void mem_create
  ( int mx, int my, int mz,
    int nx, int ny, int nz,
    int gx, int gy, int gz );

  // Groups

  /// Return the number of subgroups in the current group
  virtual int group_count () const throw();

  /// Return the name of the ith subgroup
  virtual std::string group_name (size_t i) const throw();

  /// Change group name for subsequent open or create
  virtual void group_chdir (std::string name) throw();

  /// Open the existing group from current group name
  virtual void group_open () throw();

  /// Create a new group from current group name (and open it)
  virtual void group_create () throw();

  /// Get the current group
  virtual void group_close () throw();

  /// Read a metadata item associated with the opened group
  virtual void group_read_meta
  ( void * buffer, std::string name,  int * s_type,
    int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();
  
  /// Write a metadata item associated with the opened group
  virtual void group_write_meta
  ( const void * buffer, std::string name, int type,
    int n1=1, int n2=0, int n3=0, int n4=0) throw();
  
};

#endif /* DISK_IFRIT_HPP */

