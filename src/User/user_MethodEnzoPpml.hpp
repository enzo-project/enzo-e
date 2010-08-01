// $Id: user_MethodEnzoPpml.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     user_MethodEnzoPpml.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon May 17 14:16:01 PDT 2010
/// @brief    Implementation of Enzo PPML method

#ifndef USER_METHOD_ENZO_PPML_HPP
#define USER_METHOD_ENZO_PPML_HPP

class MethodEnzoPpml : public UserMethod {

  /// @class    MethodEnzoPpml
/// @ingroup  Enzo
/// @brief    Encapsulate Enzo's PPML hydro method

public: // interface

/// Creae a new MethodEnzoPpml object
MethodEnzoPpml(Global * global,
	       Enzo * enzo)
  : UserMethod (global),
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

  void advance_block( DataBlock * data_block,
			      double t, double dt ) throw(); 

  /// Return the name of the method

  virtual std::string method_name() const throw() 
  { return "ppm"; };
  
private: // attributes

  Enzo * enzo_;

};

#endif /* USER_METHOD_ENZO_PPML_HPP */
