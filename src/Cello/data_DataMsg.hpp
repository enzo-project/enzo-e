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

  static long counter[CONFIG_NODE_SIZE];

  DataMsg() 
    : field_face_   (NULL),
      field_data_   (NULL),
      particle_data_(NULL),
      field_face_delete_   (false),
      field_data_delete_   (false),
      particle_data_delete_(false)
      
  {
    ++counter[cello::index_static()]; 
  }

  ~DataMsg()
  {
    --counter[cello::index_static()];
    
    if (field_face_delete_) {
      delete field_face_;
      field_face_ = NULL;
    }
    if (field_data_delete_) {
      delete field_data_;
      field_data_ = NULL;
    }
    if (particle_data_delete_) {
      delete particle_data_;
      particle_data_ = NULL;
    }
  }

  /// Copy constructor
  DataMsg(const DataMsg & data_msg) throw()
  {
    ++counter[cello::index_static()]; 
  };

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
  void set_field_face  (FieldFace * field_face, bool is_new) 
  {
    field_face_ = field_face; 
    field_face_delete_ = is_new;
  }

  /// Return the serialized FieldFace array
  char * field_array () 
  { return field_array_; }

  /// Return the ParticleData
  ParticleData * particle_data () 
  { return particle_data_; }

  /// Set the ParticleData object
  void set_particle_data  (ParticleData * particle_data, bool is_new) 
  {
    particle_data_ = particle_data; 
    particle_data_delete_ = is_new;
  }

  /// Delete the ParticleData object
  void delete_particle_data  () 
  { 
    delete particle_data_; particle_data_ = NULL; 
  }

  /// Return the FieldData
  FieldData * field_data () 
  { return field_data_; }

  /// Set the FieldData object
  void set_field_data    (FieldData * field_data, bool is_new) 
  {
    field_data_ = field_data;
    field_data_delete_ = is_new;
  }

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

  /// Debugging
  void print ()
  {
    CkPrintf ("%d %p TRACE_DATA_MSG field_face_    = %p\n",CkMyPe(),this,field_face_);
    CkPrintf ("%d %p TRACE_DATA_MSG field_data_    = %p\n",CkMyPe(),this,field_data_);
    CkPrintf ("%d %p TRACE_DATA_MSG particle_data_ = %p\n",CkMyPe(),this,particle_data_);
    CkPrintf ("%d %p TRACE_DATA_MSG field_face_delete_    = %d\n",CkMyPe(),this,field_face_delete_);
    CkPrintf ("%d %p TRACE_DATA_MSG field_data_delete_    = %d\n",CkMyPe(),this,field_data_delete_);
    CkPrintf ("%d %p TRACE_DATA_MSG particle_data_delete_ = %d\n",CkMyPe(),this,particle_data_delete_);
    fflush(stdout);

  }
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

  /// Whethere FieldFace data should be deleted in destructor
  bool field_face_delete_;

  /// Whethere FieldData data should be deleted in destructor
  bool field_data_delete_;

  /// Whethere FieldFace data should be deleted in destructor
  bool particle_data_delete_;

};

#endif /* DATA_DATA_MSG_HPP */

