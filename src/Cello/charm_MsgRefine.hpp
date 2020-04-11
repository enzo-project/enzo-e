// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgRefine.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Charm] Declaration of the MsgRefine Charm++ Message

#ifndef CHARM_MSG_REFINE_HPP
#define CHARM_MSG_REFINE_HPP

#include "cello.hpp"

class Block;
class Data;
class DataMsg;

class MsgRefine : public CMessage_MsgRefine {

public: // interface

  friend class Block;
  static long counter[CONFIG_NODE_SIZE];

  MsgRefine();

  MsgRefine
  (Index index,
   int nx, int ny, int nz,
   int num_field_blocks,
   int num_adapt_steps,
   int cycle, double time, double dt,
   int refresh_type,
   int num_face_level, int * face_level) ;

  virtual ~MsgRefine();

  /// Copy constructor
  MsgRefine(const MsgRefine & data_msg) throw()
  {
#ifdef DEBUG_MSG_REFINE  
    CkPrintf ("%d %s:%d DEBUG_MSG_REFINE creating %p(%p)\n",CkMyPe(),__FILE__,__LINE__,this,&data_msg);
#endif  
    ++counter[cello::index_static()]; 
  };

  /// Assignment operator
  MsgRefine & operator= (const MsgRefine & data_msg) throw()
  {
#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d %s:%d DEBUG_MSG_REFINE assigning %p = %p\n",CkMyPe(),__FILE__,__LINE__,
	    this,&data_msg);
#endif  
    return *this;
  }

  /// Set the DataMsg object
  void set_data_msg (DataMsg * data_msg);
  
  /// Update the Data with data stored in this message
  void update (Data * data);

  void print()
  {
    CkPrintf ("bool is_local_ = %d\n", is_local_);
    CkPrintf ("double time_ = %g\n", time_);
    CkPrintf ("double dt_ = %g\n", dt_);
    CkPrintf ("DataMsg * data_msg_ = %p\n", data_msg_);
    CkPrintf ("void * buffer_ = %p\n",buffer_);
    CkPrintf ("\n");
    CkPrintf ("int nx_, ny_, nz_ = %d %d %d\n",nx_, ny_, nz_);
    CkPrintf ("int num_field_blocks_ = %d\n",num_field_blocks_);
    CkPrintf ("int num_adapt_steps_ = %d\n",num_adapt_steps_);
    CkPrintf ("int cycle_ = %d\n",cycle_);
    CkPrintf ("int refresh_type_ = %d\n",refresh_type_);
    CkPrintf ("int num_face_level_ = %d\n",num_face_level_);
    CkPrintf ("int * face_level_ = %p\n",face_level_);
  }
public: // static methods

  /// Pack data to serialize
  static void * pack (MsgRefine*);

  /// Unpack data to de-serialize
  static MsgRefine * unpack(void *);
  
protected: // attributes

  /// Whether destination is local or remote
  bool is_local_;

  /// Field and particle data
  DataMsg * data_msg_;

  /// MsgRefine-specific attributes
 
  double time_;          // attribute-01
  double dt_;            // attribute-02
  Index index_;          // attribute-03
  int nx_, ny_, nz_;     // attribute-04
  int num_field_blocks_; // attribute-05
  int num_adapt_steps_;  // attribute-06
  int cycle_;            // attribute-07
  int refresh_type_;     // attribute-08
  int num_face_level_;   // attribute-09
  int * face_level_;     // attribute-10

  /// Saved Charm++ buffers for deleting after unpack()
  void * buffer_;

};

#endif /* CHARM_MSG_HPP */

