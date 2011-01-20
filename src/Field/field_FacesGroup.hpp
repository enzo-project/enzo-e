// $Id: field_FacesGroup.hpp 1688 2010-08-03 22:34:22Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FACESGROUP_HPP
#define FIELD_FACESGROUP_HPP

/// @file     field_FacesGroup.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Nov 30 16:15:31 PST 2010
/// @brief    [\ref Field] Class for storing collections of faces for communication

class FacesGroup {

  /// @class    FacesGroup
  /// @ingroup  Field
  /// @brief    [\ref Field] Class for storing collections of faces for communication

public: // interface

  /// Constructor
  FacesGroup() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~FacesGroup() throw();

  /// Copy constructor
  FacesGroup(const FacesGroup & faces_group) throw();

  /// Assignment operator
  FacesGroup & operator= (const FacesGroup & faces_group) throw();

private: // functions

private: // attributes

};

#endif /* FIELD_FACESGROUP_HPP */

