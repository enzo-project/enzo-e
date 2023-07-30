// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgOrder.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2023-07-22
/// @brief    [\ref Charm] Declaration of the MsgOrder Charm++ ABC

#ifndef CHARM_MSG_ORDER_HPP
#define CHARM_MSG_ORDER_HPP

#include "cello.hpp"

class MsgOrder : public CMessage_MsgOrder {

public: // interface

  static long counter[CONFIG_NODE_SIZE];

  MsgOrder()
    : CMessage_MsgOrder(),
      index_(0),
      count_(0),
      windex_(0.0),
      wcount_(0.0)
  {
    ++counter[cello::index_static()];
    ic3_[0] = ic3_[1] = ic3_[2] = 0;
  }

  virtual ~MsgOrder()
  {
    --counter[cello::index_static()];
  }

  /// Copy constructor
  MsgOrder(const MsgOrder & msg_order) throw()
  {
    ++counter[cello::index_static()]; 
  };

  // /// Assignment operator
  // MsgOrder & operator= (const MsgOrder & msg) throw()
  // {
  //   return *this;
  // }

  inline void get_index (int & index, double & windex) const
  {
    index  = index_;
    windex = windex_;
  }
  inline void set_index (int index, double windex)
  {
    index_ = index;
    windex_ = windex;
  }
  inline void get_count (int & count, double & wcount) const
  {
    count  = count_;
    wcount = wcount_;
  }
  inline void set_count (int count, double wcount)
  {
    count_ = count;
    wcount_ = wcount;
  }
  inline void get_child (int ic3[3]) const
  {
    ic3[0]=ic3_[0];
    ic3[1]=ic3_[1];
    ic3[2]=ic3_[2];
  }
  inline void set_child (int ic3[3])
  {
    ic3_[0]=ic3[0];
    ic3_[1]=ic3[1];
    ic3_[2]=ic3[2];
  }

public: // static methods

protected: // attributes

  int index_;
  int count_;
  double windex_;
  double wcount_;
  int ic3_[3];
};

#endif /* CHARM_MSG_HPP */

