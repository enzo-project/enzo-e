// See LICENSE_CELLO file for license and copyright information

#ifndef IO_IO_HPP
#define IO_IO_HPP

/// @file     io_Io.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Aug 25 14:32:53 PDT 2011
/// @brief    [\ref Io] Declaration of the Io class

/// @enum  file_type
/// @brief Type of file to be read or written
enum file_type {
  file_type_hdf5,
  file_type_ifrit,
  file_type_png,
  file_type_adios
};

class Io {

  /// @class    Io
  /// @ingroup  Io
  /// @brief [\ref Io] Class for writing objects to disk,
  /// acting as an interface between Simulation/Hierarchy/Patch/Block
  /// objects and Disk File objects

public: // interface

  /// Constructor
  Io() throw()
  : file_(0),
    path_(""),
    name_("") {};

  /// Destructor
  virtual ~Io() throw()
  {};

  /// Create a File of the given type
  virtual void create (file_type type) throw() = 0;

  /// Open a file
  virtual void open (const char * name, 
		     const char * mode) throw() = 0;

  /// Close a file
  virtual void close () throw() = 0;

  /// Write header for the associated object
  virtual void write_header() throw () = 0;

  /// Write data for the associated object
  virtual void write_data() throw () = 0;

protected: // functions

protected: // attributes

  File * file_;
  std::string path_;
  std::string name_;

private:   // functions

  /// Prohibit copying 
  Io(const Io & io) throw();

  /// Prohibit assignment
  Io & operator= (const Io & io) throw();


};

#endif /* IO_IO_HPP */

