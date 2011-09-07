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

  /// Constructor when creating a new File
  IoField(std::string           directory, 
	  std::string           filename, 
	  enum file_format_type format) throw();

  /// Constructor when accessing an existing File
  IoField(File & file) throw();

  /// Destructor
  virtual ~IoField() throw();

  /// Create a new file 
  virtual void create () throw();

  /// Open an existing file
  virtual void open () throw();

  /// Close the file
  virtual void close () throw();

  /// Write the header for the associated Field
  virtual void write_header() throw ();

  /// Write the data for the associated Field
  virtual void write_data() throw ();

  /// Read the header for the associated Field
  virtual void read_header() throw ();

  /// Read the data for the associated Field
  virtual void read_data() throw ();

protected: // functions

protected: // attributes

private:   // functions

  /// Prohibit default constructor
  IoField() throw();

  /// Prohibit copying 
  IoField(const IoField & io_field) throw();

  /// Prohibit assignment
  IoField & operator= (const IoField & io_field) throw();


};

#endif /* IO_IO_FIELD_HPP */

