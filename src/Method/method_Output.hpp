// $Id: method_Output.hpp 1896 2010-12-03 23:54:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Output.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Method] Declaration for the Output component

#ifndef METHOD_OUTPUT_HPP
#define METHOD_OUTPUT_HPP

class Output {

  /// @class    Output
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate functions for simulation output

public: // interface

  /// Create a new Output
  Output() throw()
  {};

public: // virtual functions

  /// Accumulate contributions of a block to field images
  virtual void images_accum (Block * block) throw();

  /// Write the field images to disk
  virtual void images_write () throw();

protected: // attributes

  
};

#endif /* METHOD_OUTPUT_HPP */
