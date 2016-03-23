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

class DataMsg {

public: // interface

  DataMsg() 
    : field_face_(NULL),
      field_data_(NULL),
      particle_data_(NULL)
  {
  }

  virtual ~DataMsg()
  {
    delete field_face_;
    field_face_ = 0;
    delete particle_data_;
    particle_data_ = 0;
  }

  /// Copy constructor
  DataMsg(const DataMsg & data_msg) throw()
  { };

  /// Assignment operator
  DataMsg & operator= (const DataMsg & data_msg) throw()
  { return *this; }

  void pup(PUP::er &p) {
    TRACEPUP;
    WARNING("DataMsg::pup()",
	    "DataMsg::pup() should never be called!");
  }

  /// Return the FieldFace
  FieldFace * field_face () 
  { return field_face_; }

  /// Set the FieldFace object
  void set_field_face  (FieldFace * field_face) 
  { field_face_ = field_face; 
  }

  /// Return the serialized FieldFace array
  char * field_array () 
  { return field_array_; }

  /// Set the FieldFace object
  void set_field_array  (char  * field_array) 
  { field_array_ = field_array; 
  }

  /// Return the ParticleData
  ParticleData * particle_data () 
  { return particle_data_; }

  /// Set the ParticleData object
  void set_particle_data  (ParticleData * particle_data) 
  { particle_data_ = particle_data; }

  /// Delete the ParticleData object
  void delete_particle_data  () 
  { 
    delete particle_data_; particle_data_ = NULL; 
  }

  /// Return the FieldData
  FieldData * field_data () 
  { return field_data_; }

  /// Set the FieldData object
  void set_field_data    (FieldData * field_data) 
  { field_data_ = field_data; }

  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer);

  /// Update the Data with the data stored in this DataMsg
  void update (Data * data, bool is_local);

public: // static methods

  
protected: // attributes

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

};

#endif /* DATA_DATA_MSG_HPP */

