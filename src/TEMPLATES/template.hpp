// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef TEMPLATE_HPP
#define TEMPLATE_HPP

/// @file     template.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    [\ref Groupname] Brief description of file template.hpp

class Classname {

  /// @class    Classname
  /// @ingroup  Groupname
  /// @brief    [\ref Groupname] Brief description of class Classname.

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


private: // attributes


};

#endif /* TEMPLATE_HPP */

