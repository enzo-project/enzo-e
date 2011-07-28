// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_GroupProcessMpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Oct 15 10:40:37 PDT 2009 
/// @brief    [\ref Parallel] Interface for the GroupProcessMpi class

#ifndef PARALLEL_GROUP_PROCESS_MPI_HPP
#define PARALLEL_GROUP_PROCESS_MPI_HPP

#ifdef CONFIG_USE_MPI

enum send_enum {
  send_standard,
  send_buffered,
  send_synchronous,
  send_ready };

enum data_enum {
  data_single,
  data_double };

class GroupProcessMpi : public GroupProcess {

  /// @class    GroupProcessMpi
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] MPI helper functions

public: // interface

  /// Initialize an GroupProcessMpi object
  GroupProcessMpi(int process_first     = 0,
		  int process_last_plus = -1) throw();

public: // interface (Group)

  /// Synchronize between all compute elements in the Group
  void barrier()  throw();

  /// Synchronize between two compute elements in the Group
  void sync (int rank, int tag=0) throw();

  //--------------------------------------------------

  /// Initiate sending an array
  void * send_begin
  (int rank_dest, void * buffer, int size, int tag=0) throw();

  /// Clean up after sending an array
  void send_end(void * handle) throw();

  /// Initiate receiving an array
  void * recv_begin
  (int rank_source, void * buffer, int size, int tag=0) throw();

  /// Clean up after receiving an array
  void recv_end(void * handle) throw();

  /// Complete sending or receiving an array
  void wait(void * handle) throw();

  /// Test completeness of sending or receiving an array
  bool test (void * handle) throw();

  //--------------------------------------------------

  /// Add an array to a list of arrays to send in bulk
  void bulk_send_add(int rank_dest, void * buffer, int size, int tag=0) throw();

  /// Initiate a bulk send of multiple arrays
  void * bulk_send() throw();

  /// Add an array to a list of arrays to receive in bulk
  void bulk_recv_add(int rank_source, void * buffer, int size, int tag=0) throw();

  /// Initiate a bulk receive of multiple arrays
  void * bulk_recv() throw();

  /// Complete a bulk send or receive of multiple arrays
  void bulk_wait(void * handle) throw();

  /// Create a Reduce object for this ProcessGroup
  Reduce * create_reduce () throw ();

  //--------------------------------------------------
  // GroupProcessMpi-specific functions
  //--------------------------------------------------

  /// Set whether send is standard, buffered, synchronous, or ready
  void set_type_send (send_enum type)  throw()
  { send_type_ = type; };

  /// Return whether send is standard, buffered, synchronous, or ready
  send_enum type_send () throw()
  { return send_type_; };

  /// Set whether send is blocking or non-blocking
  void set_send_blocking (bool blocking)  throw()
  { send_blocking_ = blocking; };

  /// Set whether send is blocking or non-blocking
  bool get_send_blocking () throw()
  { return send_blocking_; };

  /// Set whether recv is blocking or non-blocking
  void set_recv_blocking (bool blocking)  throw()
  { recv_blocking_ = blocking; };

  /// Set whether recv is blocking or non-blocking
  bool get_recv_blocking () throw()
  { return recv_blocking_; };

  MPI_Comm comm ()  const throw()
  { return comm_; };

private: // functions

  void call_mpi_(const char * file, 
		 int line , 
		 const char * name, 
		 int ierr) throw();

private: // attributes

  /// Communicator for the group
  MPI_Comm comm_;

  /// First process in the group
  int process_first_;

  /// Last process (plus one) in the group 
  int process_last_plus_;

  /// Whether to use standard, buffered, synchronous, or ready sends
  enum send_enum send_type_;

  /// Whether to use blocking sends
  bool send_blocking_;
  
  /// Whether to use blocking receives
  bool recv_blocking_;

};

#endif /* CONFIG_USE_MPI */

#endif /* PARALLEL_GROUP_PROCESS_MPI_HPP */

