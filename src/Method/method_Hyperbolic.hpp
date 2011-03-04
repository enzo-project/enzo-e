// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Hyperbolic.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Method] Declaration for the Hyperbolic class

#ifndef METHOD_HYPERBOLIC_HPP
#define METHOD_HYPERBOLIC_HPP

class Hyperbolic : public Method {

  /// @class    Hyperbolic
  /// @ingroup  Method
  /// @brief    [\ref Method] Interface to a hyperbolic method

public: // interface

  /// Create a new Hyperbolic
  Hyperbolic(Parameters * parameters) throw()
    : Method(parameters)
  { }

public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute_block( DataBlock * data_block,
			      double t, double dt ) throw() = 0; 

protected: // functions


protected: // attributes

};

#endif /* METHOD_HYPERBOLIC_HPP */
