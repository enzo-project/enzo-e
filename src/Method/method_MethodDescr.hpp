// $Id: method_MethodDescr.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef METHOD_METHOD_DESCR_HPP
#define METHOD_METHOD_DESCR_HPP

/// @file     method_MethodDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Declaration of the MethodDescr class

class MethodDescr {

  /// @class    MethodDescr
  /// @ingroup  Method
  /// @brief    Top-level description of numerical methods

public: // interface

  /// Constructor
  MethodDescr() throw()
  : methods_(0)
  {};

  /// Add a method
  void add_method (Method * method)
  { methods_.push_back(method); };

  /// Get ith method
  Method * method (int i)
  { return methods_.at(i); };

private: // attributes

  std::vector<Method *> methods_;

};

#endif /* METHOD_METHODDESCR_HPP */

