// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Data] Declaration of the DataMsg Charm++ Message

#ifndef DATA_DATA_MSG_HPP
#define DATA_DATA_MSG_HPP

class ParticleData;
class FieldData;
class FieldFace;

class DataMsg : public CMessage_DataMsg {

public: // interface

  static int id_count;
  DataMsg() 
    : CMessage_DataMsg(),
      is_local_(true),
      id_(-1),
      field_face_(NULL)
  {
    field_array_ = NULL;  
  }

  virtual ~DataMsg()
  {  }

  /// Copy constructor
  DataMsg(const DataMsg & data_msg) throw()
  { };

  /// Assignment operator
  DataMsg & operator= (const DataMsg & data_msg) throw()
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
  static void * pack (DataMsg*);

  /// Unpack data to de-serialize
  static DataMsg * unpack(void *);
  
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

