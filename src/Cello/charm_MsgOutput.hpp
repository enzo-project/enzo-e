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
class Index;
class MethodOutput;
class FileHdf5;

#define TAG_LEN 9

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
    is_local_      = msg_output.is_local_;
    index_send_    = msg_output.index_send_;
    block_trace_   = msg_output.block_trace_;
    method_output_ = msg_output.method_output_;
    file_          = msg_output.file_;
    data_msg_      = msg_output.data_msg_;
    strncpy(tag_,msg_output.tag_,TAG_LEN);
    block_name     = msg_output.block_name;
    buffer_        = nullptr;
  };

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
  void set_block (Block * block);

  /// Clear the Block data in this message
  void del_block();

  void print (const char * msg);

  const char * tag() { return tag_;}
  
public: // static methods

  /// Pack data to serialize
  static void * pack (MsgOutput*);

  /// Unpack data to de-serialize
  static MsgOutput * unpack(void *);

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

  char tag_[TAG_LEN];
  
public: // attributes

  /// Block meta-data
  std::string block_name;
  
};

#endif /* CHARM_MSG_OUTPUT_HPP */

