// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpml.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon May 17 14:16:01 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo PPML method

#ifndef ENZO_ENZO_METHOD_PPML_HPP
#define ENZO_ENZO_METHOD_PPML_HPP

class EnzoMethodPpml : public Hyperbolic {

/// @class    EnzoMethodPpml
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulate Enzo's PPML hydro method

public: // interface

/// Creae a new EnzoMethodPpml object
  EnzoMethodPpml(Error      * error,
		 Monitor    * monitor,
		 Parameters * parameters,
		 EnzoDescr * enzo)
    : Hyperbolic (error,monitor,parameters),
      enzo_(enzo)
  {};

  /// Perform any method-specific initialization

  void initialize(DataDescr * data_descr) throw();

  /// Perform any method-specific finalizations steps, e.g. to
  /// deallocate any dynamically-allocated memory

  void finalize(DataDescr * data_descr) throw();

  /// Initialize PPM variable that may change.  Called once per
  /// block per timestep.

  void initialize_block (DataBlock * data_block) throw();

  /// Finalize PPM after advancing a block a timestep, e.g. to deallocate
  /// any dynamically-allocated variables

  void finalize_block (DataBlock * data_block) throw();

  /// Apply the method to advance a block one timestep 

  void compute_block( DataBlock * data_block,
		      double t, double dt ) throw(); 

  /// Return the name of the method

  virtual std::string name() const throw() 
  { return "ppm"; };
  
private: // attributes

  EnzoDescr * enzo_;

};

#endif /* ENZO_ENZO_METHOD_PPML_HPP */
