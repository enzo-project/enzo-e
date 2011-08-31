// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_IO_HPP
#define MESH_IO_HPP

/// @file     mesh_IoField.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Aug 25 14:32:53 PDT 2011
/// @brief    [\ref Mesh] Declaration of the IoField class

class IoField {

  /// @class    IoField
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Class for writing Mesh objects, acting as an interface between Hierarchy/Patch/Block objects and Disk File objects

public: // interface

  /// Constructor
  IoField() throw()
  : path_(""),name_("") {};

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

#endif /* MESH_IO_HPP */

