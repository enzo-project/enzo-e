// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_MPI_HPP
#define PARALLEL_MPI_HPP

#include <mpi.h>

/// @file     mpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Oct 15 10:40:37 PDT 2009 
/// @bugs     Should use the Singleton design pattern
/// @brief    Interface for the ParallelMpi class

enum enum_send_type {
  send_type_standard,
  send_type_buffered,
  send_type_synchronous,
  send_type_ready };

enum enum_data_type {
  data_type_single,
  data_type_double };

class ParallelMpi {

  /// @class    ParallelMpi
  /// @ingroup  Parallel
  /// @brief    MPI helper functions

public: // interface

  /// Initialize an ParallelMpi object
  ParallelMpi();

  /// Delete an ParallelMpi object
  ~ParallelMpi();

  /// Initialize MPI
  void initialize(int * argc, char ***argv)
  {
    MPI_Init(argc,argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
  };

  /// Finalize MPI
  void finalize()
  {
    MPI_Finalize();
  }

  /// Get MPI size
  int size() 
  { return size_; }

  /// Get MPI rank
  int rank() 
  { return rank_; }

  /// Initiate sending an array
  void send_begin(char *         array, 
		  enum_data_type data_type,
		  int            dim,
		  int            * nd,
		  int            * n,
		  int            * ns = 0);
  /// Complete sending an array
  void send_end();
  /// Initiate receiving an array
  void recv_begin();
  /// Complete receiving an array
  void recv_end();

  /// Set blocking mode
  void set_send_blocking (bool send_blocking) 
  { send_blocking_ = send_blocking; }

  /// Get send_blocking mode
  bool get_send_blocking () 
  { return send_blocking_; }

  /// Set recv_blocking mode
  void set_recv_blocking (bool recv_blocking) 
  { recv_blocking_ = recv_blocking; }

  /// Get recv_blocking mode
  bool get_recv_blocking () 
  { return recv_blocking_; }

  /// Set send type: standard, buffered, synchronous, or ready
  void set_send_type (enum_send_type send_type) 
  { send_type_ = send_type; }

  /// Get send type: standard, buffered, synchronous, or ready
  enum_send_type get_send_type () 
  { return send_type_; }

private: // attributes

  /// MPI global size
  int size_;

  /// MPI global rank
  int rank_;

  /// Whether to use blocking sends
  bool send_blocking_;
  
  /// Whether to use blocking receives
  bool recv_blocking_;

  /// Whether to use standard, buffered, synchronous, or ready sends
  enum_send_type send_type_;


};

#endif /* PARALLEL_MPI_HPP */

