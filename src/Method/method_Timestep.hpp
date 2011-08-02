// See LICENSE_CELLO file for license and copyright information

/// @file     method_Timestep.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Method] Declaration for the Timestep component

#ifndef METHOD_TIMESTEP_HPP
#define METHOD_TIMESTEP_HPP


class Timestep {

  /// @class    Timestep
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate the computation of a timestep

public: // interface

  /// Create a new Timestep
  Timestep() throw()
  {};

public: // virtual functions

  /// Compute the timestep for the block

  virtual double compute ( const FieldDescr * field_descr,
			   Block * block ) throw() = 0; 

};

#endif /* METHOD_TIMESTEP_HPP */
