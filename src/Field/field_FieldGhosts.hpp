// $Id: field_FieldGhosts.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELDGHOSTS_HPP
#define FIELD_FIELDGHOSTS_HPP

/// @file     field_FieldGhosts.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Interface for the FieldGhosts class

class FieldGhosts {

  /// @class    FieldGhosts
  /// @ingroup  Field
  /// @brief    Class for representing and operating on ghost zones

public: // interface

  /// Constructor
  FieldGhosts() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~FieldGhosts() throw();

  /// Copy constructor
  FieldGhosts(const FieldGhosts & FieldGhosts) throw();

  /// Assignment operator
  FieldGhosts & operator= (const FieldGhosts & FieldGhosts) throw();

private: // functions


private: // attributes

};

#endif /* FIELD_FIELDGHOSTS_HPP */

