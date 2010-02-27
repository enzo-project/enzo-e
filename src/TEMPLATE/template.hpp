// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef FOO_HPP
#define FOO_HPP

/// @file     template.hpp
/// @brief    Brief description of file template.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010

class Foo {

  /// @class    Foo
  /// @brief    Brief description of class Foo.
  /// @ingroup  Template

public: // interface

  /// Initialize the Foo object
  Foo();

  /// Delete the Foo object
  ~Foo();

private: // functions

  /// Brief description of Foo::private_function_() in template.hpp
  private_function_();

private: // attributes

  /// Brief description of Foo::private_attribute_ in template.hpp
  int private_attribute_;

};

#endif /* FOO_HPP */

