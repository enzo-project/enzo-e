// See LICENSE_CELLO file for license and copyright information

#ifndef IO_TYPE_IMAGE_HPP
#define IO_TYPE_IMAGE_HPP

/// @file     io_TypeImage.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-02
/// @brief    [\ref Io] Declaration of the TypeImage class
///

class TypeImage {

  /// @class    TypeImage
  /// @ingroup  Io
  /// @brief [\ref Io] Class for handling the image-specific requirements
  /// of file IO.

public: // interface

  /// Constructor
  TypeImage() throw();

  /// Initialize this file Type

  virtual void initialize(Simulation * simulation,
			  File * file);

  /// Finalize this file Type

  virtual void finalize();

  /// Read the specified block from a File

  virtual void read(Hierarchy * hierarchy,
		    Patch     * patch,
		    Block     * block,
		    int         field_id);
  
  /// Write the specified block to a File

  virtual void write(const Hierarchy * hierarchy,
		     const Patch     * patch,
		     const Block     * block,
		     int               field_id);

private: // functions

    

private: // attributes


  double * values_;
};

#endif /* IO_TYPE_IMAGE_HPP */

