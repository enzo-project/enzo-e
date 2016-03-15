// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsgCoarsen.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Data] Declaration of the DataMsgCoarsen Charm++ Message

#ifndef DATA_DATA_MSG_COARSEN_HPP
#define DATA_DATA_MSG_COARSEN_HPP

class ParticleData;
class FieldData;
class FieldFace;

class DataMsgCoarsen : public CMessage_DataMsgCoarsen {

public: // interface

  static long counter;

  static int id_count;
  DataMsgCoarsen() 
    : CMessage_DataMsgCoarsen(),
      is_local_(true),
      id_(-1),
      field_face_(NULL),
      field_data_(NULL),
      particle_data_(NULL),
      buffer_(NULL)
  {  ++counter; }

  virtual ~DataMsgCoarsen();

  /// Copy constructor
  DataMsgCoarsen(const DataMsgCoarsen & data_msg) throw()
  { ++counter; };

  /// Assignment operator
  DataMsgCoarsen & operator= (const DataMsgCoarsen & data_msg) throw()
  { return *this; }

  /// Return the FieldFace
  FieldFace * field_face () 
  { return field_face_; }

  /// Set the FieldData object
  void set_field_data  (FieldData * field_data) 
  { field_data_ = field_data; }

  /// Set the ParticleData object
  void set_particle_data  (ParticleData * particle_data) 
  { particle_data_ = particle_data; }

  /// Set the FieldFace object
  void set_field_face    (FieldFace * field_face) 
  { field_face_ = field_face; }

  /// Update the Data with data stored in this message
  void update (Data * data);

public: // static methods

  /// Pack data to serialize
  static void * pack (DataMsgCoarsen*);

  /// Unpack data to de-serialize
  static DataMsgCoarsen * unpack(void *);
  
protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// Id identifying message (TEMPORARY FOR DEBUGGING)
  int id_;

  /// Field Face Data
  FieldFace * field_face_;

  /// Field data
  union {

    /// Field data if local
    FieldData * field_data_;

    /// packed source field data if remote
    char * field_array_;

  };
  
  /// Particle data
  ParticleData * particle_data_;

  /// Saved Charm++ buffer for deleting after unpack()
  void * buffer_;

};

#endif /* DATA_DATA_MSG_HPP */

