// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgRefresh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Charm] Declaration of the MsgRefresh Charm++ Message

#ifndef CHARM_MSG_REFRESH_HPP
#define CHARM_MSG_REFRESH_HPP

// #define TRACE_MSG_REFRESH /* also uncomment in control_refresh.cpp */

#include "cello.hpp"

class Data;
class DataMsg;

class MsgRefresh : public CMessage_MsgRefresh {

public: // interface

  static long counter[CONFIG_NODE_SIZE];

  MsgRefresh() ;

  virtual ~MsgRefresh();

  /// Copy constructor
  MsgRefresh(const MsgRefresh & msg_refresh) throw()
  {
    ++counter[cello::index_static()]; 
  };

  /// Assignment operator
  MsgRefresh & operator= (const MsgRefresh & data_msg) throw()
  {
    return *this;
  }

  // Set the new refresh object id
  void set_refresh_id(int id_refresh)
  { id_refresh_ = id_refresh; }

  int id_refresh() const
  { return id_refresh_; }
  
  // Set the DataMsg object
  void set_data_msg (DataMsg * data_msg);

  /// Update the Data with data stored in this message
  void update (Data * data);

  void set_block_type(std::string block, std::string type)
  {
#ifdef TRACE_MSG_REFRESH      
    name_block_ = block;
    name_type_ = type;;
#endif    
  }

#ifdef TRACE_MSG_REFRESH      
  std::string block_name ()
  {
    return name_block_;
  }
  std::string type_name()
  {
    return name_type_;
  }
#endif    

  void print(const char * message);
  
public: // static methods

  /// Pack data to serialize
  static void * pack (MsgRefresh*);

  /// Unpack data to de-serialize
  static MsgRefresh * unpack(void *);
  
protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// New Refresh object id associated with the message
  int id_refresh_;

  DataMsg * data_msg_;

  /// Saved Charm++ buffer for deleting after unpack()
  void * buffer_;

#ifdef TRACE_MSG_REFRESH      
  std::string name_block_;
  std::string name_type_;
#endif  
};

#endif /* CHARM_MSG_HPP */

