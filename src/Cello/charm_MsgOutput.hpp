// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgOutput.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-04-20
/// @brief    [\ref Charm] Declaration of the MsgOutput Charm++ Message

#ifndef CHARM_MSG_OUTPUT_HPP
#define CHARM_MSG_OUTPUT_HPP

#include "cello.hpp"
#include "mesh_BlockTrace.hpp"
class Block;
class BlockTrace;
class Data;
class DataMsg;
class FileHdf5;
class Index;
class IoBlock;
class Factory;
class MethodOutput;

class MsgOutput : public CMessage_MsgOutput {

public: // interface

  friend class Block;
  static long counter[CONFIG_NODE_SIZE];

  MsgOutput();

  MsgOutput (BlockTrace block_trace,
             MethodOutput * method_output,
             FileHdf5 * file) ;

  virtual ~MsgOutput();

  /// Copy constructor
  MsgOutput(const MsgOutput & msg_output) throw()
    : CMessage_MsgOutput() // do NOT call copy constructor on base
  {
    ++counter[cello::index_static()];
    copy_(msg_output);
    cello::hex_string(tag_,TAG_LEN); // add new tag for new message
 };

  MsgOutput & operator = (const MsgOutput & msg_output)
  {
    copy_(msg_output);
    return *this;
  }

  /// Set the DataMsg object
  void set_data_msg (DataMsg * data_msg);

  /// Copy data from this message into the provided Data object
  void update (Data * data);

  /// Return pointer to the BlockTrace object
  BlockTrace * block_trace() { return &block_trace_; }

  /// Return the Index of the sending Block
  Index index_send();

  /// Return the FileHdf5 pointer (ONLY USABLE ON WRITER)
  FileHdf5 * file() { return file_; }

  /// Set the sending index
  void set_index_send (Index index);

  /// Set the Block whose data this message will be holding
  void set_block (Block * block, const Factory * factory);

  /// Clear the Block data in this message
  void del_block();

  void print (const char * msg);

  const char * tag() { return tag_;}

  IoBlock * io_block() { return io_block_; }

  std::string block_name() const  { return block_name_; }

  double * block_lower() { return block_lower_; }
  double * block_upper() { return block_upper_; }

public: // static methods

  /// Pack data to serialize
  static void * pack (MsgOutput*);

  /// Unpack data to de-serialize
  static MsgOutput * unpack(void *);

protected: // methods

  void copy_(const MsgOutput & msg_output)
  {
    is_local_      = msg_output.is_local_;
    index_send_    = msg_output.index_send_;
    block_trace_   = msg_output.block_trace_;
    method_output_ = msg_output.method_output_;
    file_          = msg_output.file_;
    data_msg_      = msg_output.data_msg_;
    buffer_        = nullptr;
    io_block_      = msg_output.io_block_;
    block_name_     = msg_output.block_name_;
    for (int i=0; i<3; i++) {
      block_lower_[i]    = msg_output.block_lower_[i];
      block_upper_[i]    = msg_output.block_upper_[i];
    }
    strncpy(tag_,msg_output.tag_,TAG_LEN);
  }

protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// Index of the sending Block
  Index index_send_;

  /// Trace to currently active Block
  BlockTrace block_trace_;

  /// Saved pointer to MethodOutput object (on writing Block only!)
  MethodOutput * method_output_;

  /// Saved pointer to FileHdf5 object (on writing Block only!)
  FileHdf5 * file_;

  /// Associated Block data to output if any
  DataMsg * data_msg_;

  /// Saved Charm++ buffers for deleting after unpack()
  void * buffer_;

  /// Random hex tag for tracking messages for debugging
  char tag_[TAG_LEN+1];

  /// Data for the Block
  IoBlock * io_block_;

  /// Block meta-data
  std::string block_name_;

  /// Block extents
  double block_lower_[3];
  double block_upper_[3];

};

#endif /* CHARM_MSG_OUTPUT_HPP */

