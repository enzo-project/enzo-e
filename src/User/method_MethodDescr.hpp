// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef METHOD_METHOD_DESCR_HPP
#define METHOD_METHOD_DESCR_HPP

/// @file     method_MethodDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @todo     consolidate different method types into one Method* list
/// @todo     add accessor functions
/// @brief    Declaration of the MethodDescr class

class MethodDescr {

  /// @class    MethodDescr
  /// @ingroup  Method
  /// @brief    Top-level container for user-implemented numerical components

public: // interface

  /// Constructor
  MethodDescr(Global * global) throw()
    : method_control_(0),
      method_timestep_(0),
      method_method_(0),
      global_ (global)
  {

    // Set "default" method control and timestepping routines

  };

  /// Destructor
  ~MethodDescr() throw()
  {
    delete method_control_;
    delete method_timestep_;
    for (size_t i=0; i<method_method_.size(); i++) {
      delete method_method_[i];
    }
  }


protected: // functions

  /// APPLICATION INHERITENCE OVERRIDE: Create named control method.
  virtual MethodControl * create_method_control_ (std::string name_method_control) = 0;

  /// APPLICATION INHERITENCE OVERRIDE: Create named timestep method.
  virtual MethodTimestep * create_method_timestep_ (std::string name_method_timestep) = 0;

  /// APPLICATION INHERITENCE OVERRIDE: Create named hyperbolic method.
  virtual MethodMethod * create_method_method_ (std::string name_method_method) = 0;

protected: // attributes

  MethodControl *             method_control_;
  MethodTimestep *            method_timestep_;
  std::vector<MethodMethod *> method_method_;
  Global     *              global_;
};

#endif /* METHOD_METHOD_DESCR_HPP */

