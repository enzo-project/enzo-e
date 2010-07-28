// $Id: parallel_GroupProcess.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_GROUP_PROCESS_HPP
#define PARALLEL_GROUP_PROCESS_HPP

/// @file     parallel_GroupProcess.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul 22 12:36:38 PDT 2010
/// @brief    Declaration of the GroupProcess class

class GroupProcess : public Group {

  /// @class    GroupProcess
  /// @ingroup  Parallel
  /// @brief    Group of distributed memory processes

public: // interface

  /// Initialize the GroupProcess object
  GroupProcess(int size = 1, int rank = 0) throw()
    : Group(size,rank)
  {  }

  /// Initiate sending an array
  virtual void * send(int rank_dest, void * buffer, int size, int tag=0) throw() = 0;

  /// Complete sending an array
  virtual void send_wait(void * handle) throw() = 0;

  /// Initiate receiving an array
  virtual void * recv(int rank_source, void * buffer, int size, int tag=0) throw() = 0;

  /// Complete receiving an array
  virtual void recv_wait(void * handle) throw() = 0;

  /// Add an array to a list of arrays to send in bulk
  virtual void bulk_send_add(int rank_dest, void * buffer, int size, int tag=0) throw() = 0;

  /// Initiate a bulk send of multiple arrays
  virtual void * bulk_send() throw() = 0;

  /// Complete a bulk send of multiple arrays
  virtual void bulk_send_wait(void * handle) throw() = 0;

  /// Add an array to a list of arrays to receive in bulk
  virtual void bulk_recv_add(int rank_source, void * buffer, int size, int tag=0) throw() = 0;

  /// Initiate a bulk receive of multiple arrays
  virtual void * bulk_recv() throw() = 0;

  /// Complete a bulk receive of multiple arrays
  virtual void bulk_recv_wait(void * handle) throw() = 0;


private: // attributes

};

#endif /* PARALLEL_GROUP_PROCESS_HPP */

