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

class MsgOutput : public CMessage_MsgOutput {

public: // interface

  friend class Block;
  static long counter[CONFIG_NODE_SIZE];

  MsgOutput();

  MsgOutput (BlockTrace block_trace, MethodOutput * method_output) ;

  virtual ~MsgOutput();

  /// Copy constructor
  MsgOutput(const MsgOutput & data_msg) throw()
  {  ++counter[cello::index_static()]; };

  /// Set the DataMsg object
  void set_data_msg (DataMsg * data_msg);
  
  /// Copy data from this message into the provided Data object
  void update (Data * data);

  /// Return pointer to the BlockTrace object
  BlockTrace * block_trace() { return &block_trace_; }

  /// Return the Index of the sending Block
  Index index_send();

  /// Set the sending index
  void set_index_send (Index index);
  
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

  /// Associated data to output if any
  DataMsg * data_msg_;

  /// Saved Charm++ buffers for deleting after unpack()
  void * buffer_;

};

#endif /* CHARM_MSG_OUTPUT_HPP */

