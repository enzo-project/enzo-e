// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_GROUP_PROCESS_HPP
#define PARALLEL_GROUP_PROCESS_HPP

/// @file     parallel_GroupProcess.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul 22 12:36:38 PDT 2010
/// @brief    [\ref Parallel] Declaration of the GroupProcess class

class Reduce;

class GroupProcess : public ParallelGroup {

  /// @class    GroupProcess
  /// @ingroup  Parallel  
  /// @todo     Support more flexible process subsets
  /// @brief    [\ref Parallel] ParallelGroup of distributed memory processes

public: // static interface

  static GroupProcess * create (int process_first     = 0,
				int process_last_plus = -1) throw();

protected: // interface

/// Protected since GroupProcess objects must be created with create()
GroupProcess(int size = 1, int rank = 0) throw()
  : ParallelGroup(size,rank)
  {  }


public: // interface

  //--------------------------------------------------

  /// Initiate sending an array
  virtual void * send_begin
  (int rank_dest, void * buffer, int size, int tag=0) throw() = 0;

  /// Clean up after sending an array
  virtual void send_end(void * handle) throw() = 0;

  /// Initiate receiving an array
  virtual void * recv_begin
  (int rank_source, void * buffer, int size, int tag=0) throw() = 0;

  /// Clean up after receiving an array
  virtual void recv_end(void * handle) throw() = 0;

  /// Wait until the send or receive operation completes
  virtual void wait(void * handle) throw() = 0;

  /// Test completeness of sending or receiving an array
  virtual bool test (void * handle) throw() = 0;

  //--------------------------------------------------

  /// Add an array to a list of arrays to send in bulk
  virtual void bulk_send_add(int rank_dest, void * buffer, int size, int tag=0) throw() = 0;

  /// Initiate a bulk send of multiple arrays
  virtual void * bulk_send() throw() = 0;

  /// Add an array to a list of arrays to receive in bulk
  virtual void bulk_recv_add(int rank_source, void * buffer, int size, int tag=0) throw() = 0;

  /// Initiate a bulk receive of multiple arrays
  virtual void * bulk_recv() throw() = 0;

  /// Complete a bulk send or receive of multiple arrays
  virtual void bulk_wait(void * handle) throw() = 0;

  //--------------------------------------------------

  /// Create a Reduce object for this ProcessGroup
  virtual Reduce * create_reduce () throw ()= 0;

};

#endif /* PARALLEL_GROUP_PROCESS_HPP */

