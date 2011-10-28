// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoTimestepPpml.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Enzo] Declaration for the EnzoTimestepPpml component

#ifndef ENZO_ENZO_TIMESTEP_HPP
#define ENZO_ENZO_TIMESTEP_HPP

class EnzoTimestepPpml : public TimestepPpml {

  /// @class    EnzoTimestepPpml
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate determination of timestep

public: // interface

  /// Create a new EnzoTimestepPpml
  EnzoTimestepPpml() throw();

public: // virtual functions

  /// Compute the timestep for the block

  virtual double compute ( const FieldDescr * field_descr,
			   Block * block ) throw(); 

};

#endif /* ENZO_ENZO_TIMESTEP_HPP */
