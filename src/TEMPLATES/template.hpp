// See LICENSE_CELLO file for license and copyright information

#ifndef GROUPNAME_CLASSNAME_HPP
#define GROUPNAME_CLASSNAME_HPP

/// @file     groupname_Classname.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Groupname] Declaration of the Classname class
///

class Classname {

  /// @class    Classname
  /// @ingroup  Groupname
  /// @brief    [\ref Groupname] 

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

#endif /* GROUPNAME_CLASSNAME_HPP */

