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
  : sequence_(0)
  {};

  /// Add a method
  void add_method (std::string method_name)
  { Method * method = create_method_(method_name);
    if (method) sequence_.push_back(method); 
  };

  /// Get ith method
  Method * operator [] (int i)
  { return sequence_.at(i); };

private: // functions

  /// Create named method.  IMPLEMENTATION IN USER SPACE
  Method * create_method_ (std::string method_name);

private: // attributes

  std::vector<Method *> sequence_;

};

#endif /* METHOD_METHODDESCR_HPP */

