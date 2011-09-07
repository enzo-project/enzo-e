// See LICENSE_CELLO file for license and copyright information

#ifndef IO_IO_HPP
#define IO_IO_HPP

/// @file     io_Io.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Aug 25 14:32:53 PDT 2011
/// @brief    [\ref Io] Declaration of the Io class

/// @enum  file_type
/// @brief Type of file to be read or written

enum file_format_type {
  file_format_hdf5,
  file_format_ifrit,
  file_format_png,
  file_format_adios
};

class Io {

  /// @class    Io
  /// @ingroup  Io
  /// @brief [\ref Io] Class for writing objects to disk,
  /// acting as an interface between Simulation/Hierarchy/Patch/Block
  /// objects and Disk File objects

public: // interface

  /// Constructor when creating a new File
  Io(std::string directory, 
     std::string filename, 
     enum file_format_type file_format) throw();

  /// Constructor when accessing an existing File
  Io(File & file) throw();

  /// Destructor
  virtual ~Io() throw();

  /// Create a new file
  virtual void create () throw() = 0;

  /// Open an existing file
  virtual void open () throw() = 0;

  /// Close the file
  virtual void close () throw() = 0;

  /// Write the header for the associated object
  virtual void write_header() throw () = 0;

  /// Write the data for the associated object
  virtual void write_data() throw () = 0;

  /// Read the header for the associated object
  virtual void read_header() throw () = 0;

  /// Read the data for the associated object
  virtual void read_data() throw () = 0;

protected: // functions

  /// Prohibit default constructor
  Io() throw();

  /// Prohibit copying 
  Io(const Io & io) throw();

  /// Prohibit assignment
  Io & operator= (const Io & io) throw() ;

protected: // attributes

  /// File object
  File * file_;

  /// Whether the file is new or pre-existing
  bool is_file_new_;

  //  Iter * iter_;

  //  Type * type_;

};

#endif /* IO_IO_HPP */

