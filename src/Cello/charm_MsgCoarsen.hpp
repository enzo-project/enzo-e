// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgCoarsen.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Charm] Declaration of the MsgCoarsen Charm++ Message

#ifndef CHARM_MSG_COARSEN_HPP
#define CHARM_MSG_COARSEN_HPP

class ParticleData;
class FieldData;
class FieldFace;
class Data;
class DataMsg;

class MsgCoarsen : public CMessage_MsgCoarsen {

public: // interface

  static long counter;

  static int id_count;
  MsgCoarsen() ;

  virtual ~MsgCoarsen();

  /// Copy constructor
  MsgCoarsen(const MsgCoarsen & data_msg) throw()
  { ++counter; };

  /// Assignment operator
  MsgCoarsen & operator= (const MsgCoarsen & data_msg) throw()
  { return *this; }

  /// Return the FieldFace
  FieldFace * field_face ();

  /// Set the DataMsg object
  void set_data_msg  (DataMsg * data_msg);

  // /// Set the FieldData object
  // void set_field_data  (FieldData * field_data);

  // /// Set the ParticleData object
  // void set_particle_data  (ParticleData * particle_data) ;

  // /// Set the FieldFace object
  // void set_field_face    (FieldFace * field_face) ;

  /// Update the Data with data stored in this message
  void update (Data * data);

public: // static methods

  /// Pack data to serialize
  static void * pack (MsgCoarsen*);

  /// Unpack data to de-serialize
  static MsgCoarsen * unpack(void *);
  
protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// Id identifying message (TEMPORARY FOR DEBUGGING)
  int id_;

  DataMsg * data_msg_;

  /// Saved Charm++ buffer for deleting after unpack()
  void * buffer_;

};

#endif /* CHARM_MSG_HPP */

