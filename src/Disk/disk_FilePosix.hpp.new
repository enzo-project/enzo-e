// // See LICENSE_CELLO file for license and copyright information

// #ifndef DISK_FILE_POSIX_HPP
// #define DISK_FILE_POSIX_HPP

// /// @file     disk_FilePosix.hpp
// /// @author   James Bordner (jobordner@ucsd.edu)
// /// @date     2011-09-15
// /// @brief    [\ref Disk] Interface for the FilePosix class

// class FilePosix : public File {

//   /// @class    FilePosix
//   /// @ingroup  Disk
//   /// @brief    [\ref Disk] Class for writing and reading basic POSIX files

// public: // interface

//   /// Create a file with the given path and filename

//   FilePosix (std::string path, std::string name) throw();


//   // Files

//   /// Open an existing file
//   virtual void file_open () throw();

//   /// Create a new file
//   virtual void file_create () throw();

//   /// Close the file
//   virtual void file_close () throw();
  
//   /// Read a metadata item associated with the file
//   virtual void file_read_meta
//   ( void * buffer, std::string name,  enum scalar_type * s_type,
//     int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();
  
//   /// Write a metadata item associated with the file
//   virtual void file_write_meta
//   ( const void * buffer, std::string name, enum scalar_type type,
//     int n0=1, int n1=0, int n2=0, int n3=0, int n4=0) throw();
  

//   // Datasets

//   /// Open an existing dataset for reading
//   virtual void data_open
//   ( std::string name,  enum scalar_type * type,
//     int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

//   /// Create a new dataset for writing (and open it)
//   virtual void data_create
//   ( std::string name,  enum scalar_type type,
//     int n0=1, int n1=1, int n2=1, int n3=1, int n4=1) throw();

//   /// Read from the opened dataset
//   virtual void data_read (void * buffer) throw();

//   /// Write to the opened dataset
//   virtual void data_write (const void * buffer) throw();

//   /// Close the opened dataset
//   virtual void data_close () throw();


// private: // functions


//   /// Return the number of bytes in the given data type
//   size_t scalar_size_ (enum scalar_type type) const;

// private: // attributes

//   /// File descriptor
//   FILE * fp_;

//   /// Whether file is open or closed
//   bool  is_file_open_;

//   /// Type of data in the HDF5 datatype
//   scalar_type data_type_;

//   /// Dataset rank, 0 to 5
//   int data_rank_;

//   /// Dataset size
//   hsize_t data_size_[5];

//   /// Whether a dataset is open or closed
//   bool  is_data_open_;

// };

// #endif /* DISK_FILE_HDF5_HPP */

