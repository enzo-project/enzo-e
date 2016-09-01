// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgRefine.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Charm] Declaration of the MsgRefine Charm++ Message

#ifndef CHARM_MSG_REFINE_HPP
#define CHARM_MSG_REFINE_HPP

#include "cello.hpp"

class Block;
class ParticleData;
class FieldData;
class FieldFace;
class Data;
class DataMsg;

class MsgRefine : public CMessage_MsgRefine {

public: // interface

  friend class Block;
  static long counter[CONFIG_NODE_SIZE];

  MsgRefine();

  MsgRefine
  (Index index,
   int nx, int ny, int nz,
   int num_field_blocks,
   int num_adapt_steps,
   int cycle, double time, double dt,
   int refresh_type,
   int num_face_level, int * face_level) ;

  virtual ~MsgRefine();

  /// Copy constructor
  MsgRefine(const MsgRefine & data_msg) throw()
  {
    ++counter[cello::index_static()]; 
  };

  /// Assignment operator
  MsgRefine & operator= (const MsgRefine & data_msg) throw()
  { return *this; }

  /// Set the DataMsg object
  void set_data_msg (DataMsg * data_msg);
  
  /// Update the Data with data stored in this message
  void update (Data * data);

public: // static methods

  /// Pack data to serialize
  static void * pack (MsgRefine*);

  /// Unpack data to de-serialize
  static MsgRefine * unpack(void *);
  
protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// Field and particle data
  DataMsg * data_msg_;

  /// Saved Charm++ buffer for deleting after unpack()
  void * buffer_;

  /// MsgRefine-specific attributes

  Index index_;
  int nx_, ny_, nz_;
  int num_field_blocks_;
  int num_adapt_steps_;
  int cycle_;
  double time_;
  double dt_;
  int refresh_type_;
  int num_face_level_;
  int * face_level_;

};

#endif /* CHARM_MSG_HPP */

