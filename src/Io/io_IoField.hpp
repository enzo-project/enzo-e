// See LICENSE_CELLO file for license and copyright information

#ifndef IO_IO_FIELD_HPP
#define IO_IO_FIELD_HPP

/// @file     io_IoField.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Aug 25 14:32:53 PDT 2011
/// @brief    [\ref Io] Declaration of the IoField class

class IoField : Io {

  /// @class    IoField
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for writing Field objects, acting as an interface between Field objects and Disk objects

public: // interface

  /// Constructor
  IoField() throw()
  : Io() {};

  /// Destructor
  virtual ~IoField() throw()
  {};

  /// Create a File of the given type
  virtual void create (file_type type) throw();

  /// Open a file
  virtual void open (const char * name, 
		     const char * mode) throw();

  /// Close a file
  virtual void close () throw();

  /// Write header for the associated object
  virtual void write_header() throw ();

  /// Write data for the associated object
  virtual void write_data() throw ();

protected: // functions

protected: // attributes

private:   // functions

  /// Prohibit copying 
  IoField(const IoField & io_field) throw();

  /// Prohibit assignment
  IoField & operator= (const IoField & io_field) throw();


};

#endif /* IO_IO_FIELD_HPP */

