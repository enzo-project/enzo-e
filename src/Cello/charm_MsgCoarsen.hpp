// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgCoarsen.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Charm] Declaration of the MsgCoarsen Charm++ Message

#ifndef CHARM_MSG_COARSEN_HPP
#define CHARM_MSG_COARSEN_HPP

#include "cello.hpp"

class ParticleData;
class FieldData;
class FieldFace;
class Data;
class DataMsg;

class MsgCoarsen : public CMessage_MsgCoarsen {

public: // interface

  static long counter[CONFIG_NODE_SIZE];

  MsgCoarsen();

  MsgCoarsen( int num_face_level, int * face_level, int ic3[3]);

  virtual ~MsgCoarsen();

  /// Copy constructor
  MsgCoarsen(const MsgCoarsen & data_msg) throw()
  {
    ++counter[cello::index_static()]; 
  };

  /// Assignment operator
  MsgCoarsen & operator= (const MsgCoarsen & data_msg) throw()
  { return *this; }

  /// Set the DataMsg object
  void set_data_msg  (DataMsg * data_msg);

  /// Update the Data with data stored in this message
  void update (Data * data);

  /// Return the ic3_ attribute
  int * ic3() { return ic3_; }

  /// Return the face_level_ attribute
  int * face_level() { return face_level_; }

public: // static methods

  /// Pack data to serialize
  static void * pack (MsgCoarsen*);

  /// Unpack data to de-serialize
  static MsgCoarsen * unpack(void *);
  
protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// Field and particle data
  DataMsg * data_msg_;

  /// Saved Charm++ buffer for deleting after unpack()
  void * buffer_;

  /// MsgRefine-specific attributes

  int num_face_level_;
  int * face_level_;
  int ic3_[3];

};

#endif /* CHARM_MSG_HPP */

