// $Id: method_MethodUserHyperbolic.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_MethodUserHyperbolic.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    Declaration for the MethodUserHyperbolic class

#ifndef METHOD_METHOD_USER_HYPERBOLIC_HPP
#define METHOD_METHOD_USER_HYPERBOLIC_HPP

class MethodUserHyperbolic : public MethodUser {

  /// @class    MethodUserHyperbolic
  /// @ingroup  Method
  /// @brief    Encapsulate hyperbolic algorithms

public:

  MethodUserHyperbolic()
    : MethodUser (){};

  ~MethodUserHyperbolic() {};
};

#endif /* METHOD_METHOD_USER_HYPERBOLIC_HPP */
