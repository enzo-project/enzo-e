// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgInitial.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-04-20
/// @brief    [\ref Charm] Declaration of the MsgInitial Charm++ Message

#ifndef CHARM_MSG_INITIAL_HPP
#define CHARM_MSG_INITIAL_HPP

#include "cello.hpp"

class Block;
class Data;
class DataMsg;
class IoBlock;
class Factory;
class MethodInitial;

class MsgInitial : public CMessage_MsgInitial {

public: // interface

  friend class Block;
  static long counter[CONFIG_NODE_SIZE];

  MsgInitial();

  virtual ~MsgInitial();

  /// Copy constructor
  MsgInitial(const MsgInitial & msg_initial) throw()
    : CMessage_MsgInitial() // do NOT call copy constructor on base
  {
    is_local_       = msg_initial.is_local_;
    buffer_         = nullptr;
    data_type_      = msg_initial.data_type_;
    data_name_      = msg_initial.data_name_;
    data_attribute_ = msg_initial.data_attribute_;
    data_precision_ = msg_initial.data_precision_;
    data_bytes_     = 0;
    data_values_    = nullptr;
    data_delete_    = true;
    count_           = msg_initial.count_;
    // new message, so new tag
    cello::hex_string(tag_,TAG_LEN);

    for (int i=0; i<4; i++) {
      n4_[i] = msg_initial.n4_[i];
      h4_[i] = msg_initial.h4_[i];
    }
    nx_ = msg_initial.nx_;
    ny_ = msg_initial.ny_;
    nz_ = msg_initial.nz_;
    IX_ = msg_initial.IX_;
    IY_ = msg_initial.IY_;
    IZ_ = msg_initial.IZ_;
  };

  /// Copy data from this message into the provided Data object
  void update (Data * data);

  void print (const char * msg);

  const char * tag() { return tag_;}

public: // static methods

  /// Pack data to serialize
  static void * pack (MsgInitial*);

  /// Unpack data to de-serialize
  static MsgInitial * unpack(void *);

public: // methods
  
  /// Set data array for a field
  void set_field_data
  (std::string field_name,
   char * data, int data_size, int data_precision);

  /// Get data array for a field
  void get_field_data
  (std::string * field_name, char ** data, int * data_precision);

  /// Set data array for a particle attribute
  void set_particle_data
  (std::string particle_name, std::string particle_attribute,
   char * data, int data_size, int data_precision);
  
  /// Set data array for a particle attribute
  void get_particle_data
  (std::string * particle_name, std::string * particle_attribute,
   char ** data, int * data_size, int * data_precision);

  /// Set dataset sizes
  void set_dataset (int n4[4], double h4[4],
                    int nx, int ny, int nz,
                    int IX, int IY, int IZ);

  /// Get dataset sizes
  void get_dataset (int n4[4], double h4[4],
                    int * nx, int * ny, int * nz,
                    int * IX, int * IY, int * IZ);
  
  /// Signal that this is the last message to be received
  void set_count (int count)
  { count_ = count; }

  int count() const { return count_; }

  std::string data_type() const { return data_type_; }
  
protected: // methods

  void copy_data_( char * data, int data_size, int data_precision);

protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// Type of data, field or particle
  std::string data_type_;

  /// Name of the field or particle type
  std::string data_name_;

  /// Particle attribute if data_type_ is "particle"
  std::string data_attribute_;

  /// Precision of data
  int data_precision_;

  /// Number of elements in data array (of given precision type)
  int data_bytes_;

  /// Data values in a packed array of length data_bytes_
  char * data_values_;

  /// Whether to delete data_values_ when deleting this message
  int data_delete_;
  
  /// If positive, the expected count of number of messages (including
  /// this one) that will be received
  int count_;

  /// Saved Charm++ buffers for deleting after unpack()
  void * buffer_;

  /// Random hex tag for tracking messages for debugging
  char tag_[TAG_LEN+1];

  /// File dataset size
  int n4_[4];
  /// Cell width (for adjusting particle displacements)
  double h4_[4];
  /// Cello data size
  int nx_,ny_,nz_;
  /// Axis remapping hdf5[4] to cello[3]
  int IX_,IY_,IZ_;

public: // attributes

};

#endif /* CHARM_MSG_INITIAL_HPP */

