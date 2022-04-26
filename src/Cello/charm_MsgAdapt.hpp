// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgAdapt.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-12-18
/// @brief    [\ref Charm] Declaration of the MsgAdapt Charm++ Message

#ifndef CHARM_MSG_ADAPT_HPP
#define CHARM_MSG_ADAPT_HPP

#include "cello.hpp"

class Adapt;
class Data;
class DataMsg;
class FieldData;
class FieldFace;
class ParticleData;

class MsgAdapt : public CMessage_MsgAdapt {

public: // interface

  friend class Block;

  static long counter[CONFIG_NODE_SIZE];

  MsgAdapt
  (
   int adapt_step,
   Index index,
   int ic3[3],
   int of3[3],
   int level_now,
   int level_min,
   int level_max,
   bool can_coarsen
   )
    :
    adapt_step_(adapt_step),
    index_(index),
    level_now_(level_now),
    level_min_(level_min),
    level_max_(level_max),
    can_coarsen_(can_coarsen),
    buffer_(nullptr)
  {
    cello::hex_string(tag_,TAG_LEN); // add new tag for new message
    ic3_[0] = ic3[0];
    ic3_[1] = ic3[1];
    ic3_[2] = ic3[2];
    ofv_[0].resize(1);
    ofv_[1].resize(1);
    ofv_[2].resize(1);
    ofv_[0][0] = of3[0];
    ofv_[1][0] = of3[1];
    ofv_[2][0] = of3[2];
    ++counter[cello::index_static()]; 
  }

  void add_face(int of3[3])
  {
    ofv_[0].push_back(of3[0]);
    ofv_[1].push_back(of3[1]);
    ofv_[2].push_back(of3[2]);
  }

    virtual ~MsgAdapt();

  /// Copy constructor
  MsgAdapt(const MsgAdapt & data_msg) throw()
  {
    ++counter[cello::index_static()]; 
  };

  /// Assignment operator
  MsgAdapt & operator= (const MsgAdapt & data_msg) throw()
  { return *this; }

  const char * tag() { return tag_;}

public: // static methods

  /// Pack data to serialize
  static void * pack (MsgAdapt*);
  
  /// Unpack data to de-serialize
  static MsgAdapt * unpack(void *);

  void print (const char * msg);
  
protected: // methods

  MsgAdapt() { }

protected: // attributes

  int adapt_step_;
  Index index_;
  int ic3_[3];
  std::vector<int> ofv_[3];
  int level_now_;
  int level_min_;
  int level_max_;
  bool can_coarsen_;
  /// Saved Charm++ buffer for deleting after unpack()
  void * buffer_;
  /// Random hex tag for tracking messages for debugging
  char tag_[TAG_LEN+1];

};

#endif /* CHARM_MSG_HPP */
