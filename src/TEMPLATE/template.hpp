// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef FILENAME_HPP
#define FILENAME_HPP

/// @file     filename.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file filename.hpp

class Classname {

  /// @class    Classname
  /// @ingroup  Groupname
  /// @brief    Brief description of class Classname.

public: // interface

  /// Constructor
  Classname() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Classname() throw();

  /// Copy constructor
  Classname(const Classname & classname) throw();

  /// Assignment operator
  Classname & operator= (const Classname & classname) throw();

private: // functions

  /// Brief description of Classname::private_function_() in filename.hpp
  private_function_() throw();

private: // attributes

  /// Brief description of Classname::private_attribute_ in filename.hpp
  int private_attribute_;

};

#endif /* FILENAME_HPP */

