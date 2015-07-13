// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_FieldMsg.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Mesh] Declaration of the FieldMsg Charm++ Message

#ifndef MESH_FIELD_MSG_HPP
#define MESH_FIELD_MSG_HPP

class FieldMsg : public CMessage_FieldMsg {

public:
  /// Array length
  int n;

  /// Array data
  char * a;

  /// Child indices
  int ic3[3];
};

#endif /* MESH_FIELD_MSG_HPP */

