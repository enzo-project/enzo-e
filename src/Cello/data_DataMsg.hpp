// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Data] Declaration of the DataMsg Charm++ Message

#ifndef DATA_DATA_MSG_HPP
#define DATA_DATA_MSG_HPP

class FluxData;
class ParticleData;
class FieldData;
class FieldFace;

class DataMsg {

public: // interface

  static long counter[CONFIG_NODE_SIZE];

  DataMsg() 
    : field_face_(nullptr),
      field_face_delete_   (false),
      field_data_(nullptr),
      field_data_delete_   (false),
      particle_data_(nullptr),
      particle_data_delete_(false),
      flux_data_(nullptr),
      flux_data_delete_(false)
  {
    ++counter[cello::index_static()]; 
  }

  ~DataMsg()
  {
    --counter[cello::index_static()];
    
    if (field_face_delete_) {
      delete field_face_;
      field_face_ = nullptr;
    }
    if (field_data_delete_) {
      delete field_data_;
      field_data_ = nullptr;
    }
    if (particle_data_delete_) {
      delete particle_data_;
      particle_data_ = nullptr;
    }
    if (flux_data_delete_) {
      delete flux_data_;
      flux_data_ = nullptr;
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
  {
#ifdef DEBUG_DATA_MSG
    CkPrintf ("%d %s:%d DEBUG_DATA_MSG %p get_field_face(%p,%d)\n",
	      CkMyPe(),__FILE__,__LINE__,field_face_,field_face_delete_);
#endif    
    return field_face_;
  }

  /// Set the FieldFace object
  void set_field_face  (FieldFace * field_face, bool is_new) 
  {
    field_face_ = field_face; 
    field_face_delete_ = is_new;
#ifdef DEBUG_DATA_MSG
    CkPrintf ("%d %s:%d DEBUG_DATA_MSG %p set_field_face(%p,%d)\n",
	      CkMyPe(),__FILE__,__LINE__,field_face_,field_face_delete_);
#endif    
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
    delete particle_data_; particle_data_ = nullptr; 
  }

  /// Return the FluxData
  FluxData * flux_data () 
  { return flux_data_; }

  /// Set the FluxData object
  void set_flux_data  (FluxData * flux_data, bool is_new) 
  {
    flux_data_ = flux_data; 
    flux_data_delete_ = is_new;
  }

  /// Delete the FluxData object
  void delete_flux_data  () 
  { 
    delete flux_data_; flux_data_ = nullptr; 
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
    const int ip=CkMyPe();
    char buf[81];
    snprintf (buf,80,"%d %p TRACE_DATA_MSG TRACE_DATA_MSG",ip,this);
    CkPrintf ("%s field_face_    = %p\n",buf,field_face_);
    CkPrintf ("%s field_data_    = %p\n",buf,field_data_);
    CkPrintf ("%s particle_data_ = %p\n",buf,particle_data_);
    CkPrintf ("%s particle_data_delete_ = %d\n",buf,particle_data_delete_);
    CkPrintf ("%s flux_data_         = %p\n",buf,flux_data_);
    CkPrintf ("%s flux_data_delete   = %d\n",buf,flux_data_delete_);
    CkPrintf ("%s field_face_delete_ = %d\n",buf,field_face_delete_);
    CkPrintf ("%s field_data_delete_ = %d\n",buf,field_data_delete_);
    fflush(stdout);

  }
public: // static methods

  
protected: // attributes

  /// Field Face Data
  FieldFace * field_face_;
  /// Whethere FieldFace data should be deleted in destructor
  bool field_face_delete_;

  /// Field data
  union {

    /// Field data if local
    FieldData * field_data_;

    /// packed source field data if remote
    char * field_array_;

  };
  /// Whethere FieldData data should be deleted in destructor
  bool field_data_delete_;
  
  /// Particle data
  ParticleData * particle_data_;
  /// Whethere Particle data should be deleted in destructor
  bool particle_data_delete_;

  /// Flux data
  FluxData * flux_data_;
  /// Whether Flux data should be deleted in destructor
  bool flux_data_delete_;

};

#endif /* DATA_DATA_MSG_HPP */

