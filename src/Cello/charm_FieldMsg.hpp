// See LICENSE_CELLO file for license and copyright information

/// @file     charm_FieldMsg.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Charm] Declaration of the FieldMsg Charm++ Message

#ifndef CHARM_FIELD_MSG_HPP
#define CHARM_FIELD_MSG_HPP

class FieldMsg : public CMessage_FieldMsg {

public: // attributes

  int child_index() const { return ic3[0] + 2*(ic3[1] + 2*(ic3[2])); }
  
  /// Array length
  int n;

  /// Array data
  char * a;

  /// Child indices
  int ic3[3];
};

#endif /* CHARM_FIELD_MSG_HPP */

