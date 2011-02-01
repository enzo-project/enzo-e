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
  Hyperbolic(Error      * error,
	     Monitor    * monitor,
	     Parameters * parameters) throw()
    : Method(error,monitor,parameters)
  { }

public: // virtual functions

  /// Perform any method-specific initialization

  virtual void initialize (DataDescr * data_descr) throw() = 0;

  /// Perform any method-specific finalizations steps, e.g. to
  /// deallocate any dynamically-allocated memory

  virtual void finalize (DataDescr * data_descr) throw() = 0;

  /// Initialize PPM variable that may change.  Called once per
  /// block per timestep.

  virtual void initialize_block (DataBlock * data_block) throw() = 0;

  /// Finalize PPM after advancing a block a timestep, e.g. to deallocate
  /// any dynamically-allocated variables

  virtual void finalize_block (DataBlock * data_block) throw() = 0;

  /// Apply the method to advance a block one timestep 

  virtual void compute_block( DataBlock * data_block,
			      double t, double dt ) throw() = 0; 

  /// Return the name of the method

  virtual std::string name() const throw() = 0;

protected: // functions

  /// Specify a field or particle type and its access type

  void add_argument_(argument_enum type,
		     std::string   name,
		     access_enum   access_type,
		     DataDescr   * data_descr = 0) throw();

protected: // attributes

};

#endif /* METHOD_HYPERBOLIC_HPP */
