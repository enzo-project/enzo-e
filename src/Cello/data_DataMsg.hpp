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

  static int id_count;
  DataMsg() 
    : is_local_(true),
      n_ff_(0),
      n_fa_(0),
      id_(-1),
      field_face_(NULL)
  { field_array_ = NULL;  }
  virtual ~DataMsg()
  {  }
  /// Copy constructor
  DataMsg(const DataMsg & data_msg) throw()
  { };

  /// Assignment operator
  DataMsg & operator= (const DataMsg & data_msg) throw()
  { return *this; }

  FieldFace * field_face () { return field_face_; }
  /// Set the FieldFace object
  void set_field_face    (FieldFace * field_face) { field_face_ = field_face; }
  /// Set the FieldData object
  void set_field_data    (FieldData * field_data) { field_data_ = field_data; }

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

  /// FieldFace array size
  int n_ff_;

  /// FieldData array size
  int n_fa_;

  /// Id identifying message (TEMPORARY FOR DEBUGGING)
  int id_;

    /// Field Face Data
    FieldFace * field_face_;

  union {

    /// Field data if local
    FieldData * field_data_;
    /// packed source field data if remote
    char * field_array_;

  };
  
  // /// Particle data
  // union {
  //   // source particle data if local
  //   ParticleData * particle_data_;
  //   // packed source particle data if remote
  //   char * particle_array_;
  // };

  /// Saved Charm++ buffer for deleting after unpack()
  void * buffer_;

};

#endif /* DATA_DATA_MSG_HPP */

