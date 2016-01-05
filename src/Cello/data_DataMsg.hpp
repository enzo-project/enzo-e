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

public:

  DataMsg() : is_local_(true), field_data_(NULL), particle_data_(NULL)
  { }

  void set_field_face    (FieldFace * field_face) { field_face_ = field_face; }

  /// Pack data to serialize
  static void * pack (DataMsg*);

  /// Unpack data to de-serialize
  static DataMsg * unpack(void *);

  /// Update the Block Data with data stored in this message, ghost or
  /// interior
  void update (Data * data);

protected:

  /// Whether destination is local or remote
  bool is_local_;

  /// Field Face Data
  FieldFace * field_face_;

  /// Field data
  union {
    // source field data if local
    FieldData * field_data_;
    // packed source field data if remote
    char * field_array_;
  };

  /// Particle data
  union {
    // source particle data if local
    ParticleData * particle_data_;
    // packed source particle data if remote
    char * particle_array_;
  };

};

#endif /* DATA_DATA_MSG_HPP */

