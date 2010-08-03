// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_GROUP_PROCESS_HPP
#define PARALLEL_GROUP_PROCESS_HPP

/// @file     parallel_GroupProcess.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul 22 12:36:38 PDT 2010
/// @brief    Declaration of the GroupProcess class

class GroupProcess : public Group {

  /// @class    Group Process
  /// @ingroup  Parallel  
  /// @todo     Support more flexible process subsets
  /// @brief    Group of distributed memory processes

public: // interface

/// Initialize the GroupProcess object
GroupProcess(int size = 1, int rank = 0) throw()
  : Group(size,rank),
    send_blocking_    (true),
    recv_blocking_    (true)
  {  }

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

  /// Set whether send is blocking or non-blocking
  void set_send_blocking (bool blocking)  throw()
  { send_blocking_ = blocking; };

  /// Set whether send is blocking or non-blocking
  bool send_blocking ()  throw()
  { return send_blocking_; };

  /// Set whether recv is blocking or non-blocking
  void set_recv_blocking (bool blocking)  throw()
  { recv_blocking_ = blocking; };

  /// Set whether recv is blocking or non-blocking
  bool recv_blocking ()  throw()
  { return recv_blocking_; };

protected: // attributes

  /// Whether to use blocking sends
  bool send_blocking_;
  
  /// Whether to use blocking receives
  bool recv_blocking_;


};

#endif /* PARALLEL_GROUP_PROCESS_HPP */

