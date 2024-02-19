// See LICENSE_CELLO file for license and copyright information

/// @file     charm_EnzoMsgCheck.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-19
/// @brief    [\ref Charm] Declaration of the EnzoMsgCheck Charm++ Message

#ifndef CHARM_ENZO_MSG_CHECK_HPP
#define CHARM_ENZO_MSG_CHECK_HPP

#include "enzo.hpp"
#include "charm_enzo.hpp"

#define ADAPT_BUFFER_SIZE 1000

class EnzoBlock;

class EnzoMsgCheck : public CMessage_EnzoMsgCheck {

public: // interface

  friend class EnzoBlock;
  friend class IoEnzoWriter;
  friend class IoEnzoReader;

  static long counter[CONFIG_NODE_SIZE];

  EnzoMsgCheck();

  virtual ~EnzoMsgCheck();

  /// Copy constructor
  EnzoMsgCheck(const EnzoMsgCheck & enzo_msg_check) throw()
    : CMessage_EnzoMsgCheck() // do NOT call copy constructor on base
  {
    ++counter[cello::index_static()];
    copy_(enzo_msg_check);
    cello::hex_string(tag_,TAG_LEN); // add new tag for new message
  }

  EnzoMsgCheck & operator = (const EnzoMsgCheck & enzo_msg_check)
  {
    copy_(enzo_msg_check);
    return *this;
  }

  /// Set the DataMsg object
  void set_data_msg (DataMsg * data_msg);

  /// Copy data from this message into the provided Data object
  void update (Data * data);

  /// Copy Block data from message into the provided Block
  void update (EnzoBlock * block);

  /// Return the Index of the sending Block
  Index index_send();

  /// Set the sending index
  void set_index_send (Index index);

  /// Set the Block whose data this message will be holding
  void set_block (Block * block);

  /// Clear the Block data in this message
  void del_block();

  void print (const char * msg);

  const char * tag() { return tag_;}

  std::string block_name() const  { return block_name_; }
  int block_level() const  { return block_level_; }

  double * block_lower() { return block_lower_; }
  double * block_upper() { return block_upper_; }
  int * block_size() { return block_size_; }

  void set_parameters (Index index_this,
                       Index index_next,
                       std::string name_this,
                       std::string name_next,
                       long long index_block,
                       bool is_first,
                       bool is_last)
  {
    index_this_ =  index_this;
    index_next_ = index_next;
    name_this_ = name_this;
    name_next_ = name_next;
    index_block_ = index_block;
    is_first_ = is_first;
    is_last_ = is_last;
  }

  void get_parameters (Index & index_this,
                       Index & index_next,
                       std::string & name_this,
                       std::string & name_next,
                       long long & index_block,
                       bool & is_first,
                       bool & is_last,
                       std::string & name_dir)
  {
    index_this =  index_this_;
    index_next = index_next_;
    name_this = name_this_;
    name_next = name_next_;
    index_block = index_block_;
    is_first = is_first_;
    is_last = is_last_;
    name_dir = name_dir_;
  }

  void set_adapt (const Adapt & adapt) {
    int size = adapt.data_size();
    ASSERT2 ("EnzoMsgCheck::set_adapt()",
             "ADAPT_BUFFER_SIZE %d is too small for Adapt object of size %d",
             sizeof(int)*ADAPT_BUFFER_SIZE,size,
             size <= sizeof(int)*ADAPT_BUFFER_SIZE);
    adapt.save_data((char *)((int*)adapt_buffer_));
  }

  void get_adapt (Adapt & adapt) {
    adapt.load_data((char *)((int*)adapt_buffer_));
  }
  void set_name_dir (std::string name_dir)
  { name_dir_ = name_dir; }

  void set_io_block(IoEnzoBlock * io_block) { io_block_ = io_block; }
  IoEnzoBlock * io_block() { return io_block_; }

public: // static methods

  /// Pack data to serialize
  static void * pack (EnzoMsgCheck*);

  /// Unpack data to de-serialize
  static EnzoMsgCheck * unpack(void *);

protected: // methods

  int size_();
  char * load_(char *);
  char * save_(char *);

  void copy_(const EnzoMsgCheck & enzo_msg_check)
  {
    is_local_      = enzo_msg_check.is_local_;
    index_send_    = enzo_msg_check.index_send_;
    data_msg_      = enzo_msg_check.data_msg_;
    buffer_        = nullptr;
    block_name_    = enzo_msg_check.block_name_;
    block_level_   = enzo_msg_check.block_level_;
    for (int i=0; i<3; i++) {
      block_lower_[i]    = enzo_msg_check.block_lower_[i];
      block_upper_[i]    = enzo_msg_check.block_upper_[i];
      block_size_[i]     = enzo_msg_check.block_size_[i];
    }
    strncpy(tag_,enzo_msg_check.tag_,TAG_LEN);
    io_block_    = enzo_msg_check.io_block_;
    index_this_  = enzo_msg_check.index_this_;
    index_next_  = enzo_msg_check.index_next_;
    name_this_   = enzo_msg_check.name_this_;
    name_next_   = enzo_msg_check.name_next_;
    index_block_ = enzo_msg_check.index_block_;
    is_first_    = enzo_msg_check.is_first_;
    is_last_     = enzo_msg_check.is_last_;
    name_dir_    = enzo_msg_check.name_dir_;
    std::copy_n ( enzo_msg_check.adapt_buffer_, ADAPT_BUFFER_SIZE,
                  adapt_buffer_);
    index_file_  = enzo_msg_check.index_file_;

    index_order_ = enzo_msg_check.index_order_;
    count_order_ = enzo_msg_check.count_order_;
  }

protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// Index of the sending Block
  Index index_send_;

  /// Associated Block data to output if any
  DataMsg * data_msg_;

  /// Saved Charm++ buffers for deleting after unpack()
  void * buffer_;

  /// Random hex tag for tracking messages for debugging
  char tag_[TAG_LEN+1];

  /// Data for the Block
  IoEnzoBlock * io_block_;

  /// Block meta-data
  std::string block_name_;

  /// Level of the block
  int block_level_;

  /// Block extents
  double block_lower_[3];
  double block_upper_[3];
  /// Block size
  int block_size_[3];

  // IoEnzoWriter::p_write() parameters
  Index index_this_;
  Index index_next_;
  std::string name_this_;
  std::string name_next_;
  long long index_block_;
  bool is_first_;
  bool is_last_;

  std::string name_dir_;

  /// Array holding serialized Array object
  int adapt_buffer_[ADAPT_BUFFER_SIZE];

  /// Index for io_reader for restart
  int index_file_;

  /// index/count for load balancing
  long long index_order_;
  long long count_order_;
};

#endif /* CHARM_ENZO_MSG_CHECK_HPP */

