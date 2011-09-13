// See LICENSE_CELLO file for license and copyright information

#ifndef IO_TYPE_HPP
#define IO_TYPE_HPP

/// @file     io_Type.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-02
/// @brief    [\ref Io] Declaration of the Type class
///

class Type {

  /// @class    Type
  /// @ingroup  Io
  /// @brief [\ref Io] Abstract base class for defining the "type" of
  /// output file (data, restart, image, etc) and perforing any
  /// pre-processing steps

public: // interface

  /// Constructor
  Type() throw();

  /// Initialize this file Type

  virtual void initialize(Simulation * simulation,
			  File      * file ) = 0;

  /// Finalize this file Type

  virtual void finalize() = 0;

  /// Read the specified block from a File

  virtual void read(Hierarchy * hierarchy,
		    Patch     * patch,
		    Block     * block,
		    int         field_id) = 0;
  
  /// Write the specified block to a File

  virtual void write(const Hierarchy * hierarchy,
		     const Patch     * patch,
		     const Block     * block,
		     const int         field_id) = 0;

private: // functions


private: // attributes


};

#endif /* IO_TYPE_HPP */

