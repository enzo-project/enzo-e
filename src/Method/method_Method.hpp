// $Id: method_Method.hpp 1942 2011-01-20 00:53:45Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Method.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Method] Declaration for the Method class

#ifndef METHOD_METHOD_HPP
#define METHOD_METHOD_HPP

class Method {

  /// @class    Method
  /// @ingroup  Method
  /// @brief    [\ref Method] Interface to an application method / analysis / visualization function.

public: // interface

  /// Create a new Method
  Method(Parameters * parameters) throw()
  {};

public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute_block( Block * block,
			      double t, double dt ) throw() = 0; 

protected: // functions


protected: // attributes

};

#endif /* METHOD_METHOD_METHOD_HPP */
