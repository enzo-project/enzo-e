// See LICENSE_CELLO file for license and copyright information

/// @file     charm_EnzoMsgInferCreate.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-10-05
/// @brief    [\ref Charm] Declaration of the EnzoMsgInferCreate Charm++ Message

#ifndef CHARM_ENZO_MSG_INFER_CREATE_HPP
#define CHARM_ENZO_MSG_INFER_CREATE_HPP

#include "enzo.hpp"
#include "charm_enzo.hpp"

#define ADAPT_BUFFER_SIZE 800

class EnzoBlock;

class EnzoMsgInferCreate : public CMessage_EnzoMsgInferCreate {

public: // interface

  static long counter[CONFIG_NODE_SIZE];

  EnzoMsgInferCreate();

  virtual ~EnzoMsgInferCreate();

  /// Copy constructor
  EnzoMsgInferCreate(const EnzoMsgInferCreate & msg) throw()
    : CMessage_EnzoMsgInferCreate() // do NOT call copy constructor on base
  {
    ++counter[cello::index_static()];
    copy_(msg);
    cello::hex_string(tag_,TAG_LEN); // add new tag for new message
  }

  EnzoMsgInferCreate & operator = (const EnzoMsgInferCreate & msg)
  {
    copy_(msg);
    return *this;
  }

  const char * tag() { return tag_;}

  void print (const char * msg);

public: // static methods

  /// Pack data to serialize
  static void * pack (EnzoMsgInferCreate*);

  /// Unpack data to de-serialize
  static EnzoMsgInferCreate * unpack(void *);

protected: // methods

  void copy_(const EnzoMsgInferCreate & msg)
  {
    is_local_ = msg.is_local_;
    buffer_        = nullptr;
    strncpy(tag_,msg.tag_,TAG_LEN);
  }

protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// Saved Charm++ buffers for deleting after unpack()
  void * buffer_;

  /// Random hex tag for tracking messages for debugging
  char tag_[TAG_LEN+1];
};

#endif /* CHARM_ENZO_MSG_INFER_CREATE_HPP */

