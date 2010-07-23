// $Id: parallel_ParallelMpi.hpp 1633 2010-07-21 18:47:53Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Mpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Oct 15 10:40:37 PDT 2009 
/// @brief    Interface for the ParallelMpi class

#ifndef PARALLEL_PARALLEL_MPI_HPP
#define PARALLEL_PARALLEL_MPI_HPP

enum send_type {
  send_standard,
  send_buffered,
  send_synchronous,
  send_ready };

enum data_type {
  data_single,
  data_double };

class GroupProcessMpi : public GroupProcess {

  /// @class    GroupProcessMpi
  /// @ingroup  Parallel
  /// @brief    MPI helper functions

public: // interface

  /// Initialize an GroupProcessMpi object (singleton design pattern)
  GroupProcessMpi() throw();

public: // interface (Group)

  /// Number of compute elements in the Group
  int size() { return size_; };

  /// Rank of the compute element in the Group
  virtual int rank() { return 0; };

  /// Synchronize between all compute elements in the Group
  virtual void barrier() { };

  /// Synchronize between two compute elements in the Group
  virtual void wait() { };

public: // interface (Group)

  /// Initiate sending an array
  int send(int rank_dest, char * buffer, int size) throw();

  /// Complete sending an array
  void send_wait(int handle) throw();

  /// Initiate receiving an array
  int recv(int rank_source, char * buffer, int size) throw();

  /// Complete receiving an array
  void recv_wait(int handle) throw();

  /// Add an array to a list of arrays to send in bulk
  void bulk_send_add(int rank_dest, char * buffer, int size) throw();

  /// Initiate a bulk send of multiple arrays
  int bulk_send() throw();

  /// Complete a bulk send of multiple arrays
  void bulk_send_wait(int handle) throw();

  /// Add an array to a list of arrays to receive in bulk
  void bulk_recv_add(int rank_source, char * buffer, int size) throw();

  /// Initiate a bulk receive of multiple arrays
  int bulk_recv() throw();

  /// Complete a bulk receive of multiple arrays
  void bulk_recv_wait(int handle) throw();

private: // attributes

  /// Whether to use blocking sends
  bool send_blocking_;
  
  /// Whether to use blocking receives
  bool recv_blocking_;

  /// Whether to use standard, buffered, synchronous, or ready sends
  enum send_type send_type_;

};

#endif /* PARALLEL_PARALLEL_MPI_HPP */

