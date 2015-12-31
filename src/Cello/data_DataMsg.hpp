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

  DataMsg() : pd_(NULL), fd_(NULL), ff_(NULL)
  { }

  void set_particle_data (ParticleData * pd) { pd_ = pd; }
  void set_field_data    (FieldData * fd) { fd_ = fd; }
  void set_field_face    (FieldFace * ff) { ff_ = ff; }

  /// Pack data to serialize
  static void * pack (DataMsg*);

  /// Unpack data to de-serialize
  static DataMsg * unpack(void *);


  /// Set Particle data
protected:

  /// Particle Data
  ParticleData * pd_;

  /// Field Data
  FieldData * fd_;

  /// Field Face Data
  FieldFace * ff_;
};

#endif /* DATA_DATA_MSG_HPP */

