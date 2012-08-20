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

typedef int send_type ;

// enum data_enum {
//   data_single,
//  data_double };

class GroupProcessMpi : public GroupProcess {

  /// @class    GroupProcessMpi
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] MPI helper functions

public: // interface

  /// Initialize an GroupProcessMpi object
  GroupProcessMpi(int process_first     = 0,
		  int process_last_plus = -1) throw();

  /// Destructor
  virtual ~GroupProcessMpi() throw()
  {}

public: // interface (Group)

  /// Synchronize between all compute elements in the Group
  void barrier()  const throw();

  /// Synchronize between two compute elements in the Group
  void sync (int rank, int tag=0) const throw();

  //--------------------------------------------------

  /// Initiate sending an array
  void * send_begin
  (int rank_dest, void * buffer, int size, int tag=0) const throw();

  /// Clean up after sending an array
  void send_end(void * handle) const throw();

  /// Initiate receiving an array
  void * recv_begin
  (int rank_source, void * buffer, int size, int tag=0) const throw();

  /// Clean up after receiving an array
  void recv_end(void * handle) const throw();

  /// Send and receive data between two processes
  void send_recv
  (int rank, void * buffer, int size, int tag=0) const throw();

  /// Complete sending or receiving an array
  void wait(void * handle) const throw();

  /// Test completeness of sending or receiving an array
  bool test (void * handle) const throw();

  //--------------------------------------------------

  /// Add an array to a list of arrays to send in bulk
  void bulk_send_add
  (int rank_dest, void * buffer, int size, int tag=0) const throw();

  /// Initiate a bulk send of multiple arrays
  void * bulk_send() const throw();

  /// Add an array to a list of arrays to receive in bulk
  void bulk_recv_add
  (int rank_source, void * buffer, int size, int tag=0) const throw();

  /// Initiate a bulk receive of multiple arrays
  void * bulk_recv() const throw();

  /// Complete a bulk send or receive of multiple arrays
  void bulk_wait(void * handle) const throw();

  /// Create a Reduce object for this ProcessGroup
  Reduce * create_reduce () const throw ();

  //--------------------------------------------------
  // GroupProcessMpi-specific functions
  //--------------------------------------------------

  /// Set whether send is standard, buffered, synchronous, or ready
  void set_type_send (send_type type)  throw()
  { send_type_ = type; };

  /// Return whether send is standard, buffered, synchronous, or ready
  send_type type_send () const throw()
  { return send_type_; };

  /// Set whether send is blocking or non-blocking
  void set_send_blocking (bool blocking)  throw()
  { send_blocking_ = blocking; };

  /// Set whether send is blocking or non-blocking
  bool get_send_blocking () const throw()
  { return send_blocking_; };

  /// Set whether recv is blocking or non-blocking
  void set_recv_blocking (bool blocking)  throw()
  { recv_blocking_ = blocking; };

  /// Set whether recv is blocking or non-blocking
  bool get_recv_blocking () const throw()
  { return recv_blocking_; };

  MPI_Comm comm () const throw()
  { return comm_; };

private: // functions

  void call_mpi_(const char * file, 
		 int line , 
		 const char * name, 
		 int ierr) const throw();

private: // attributes

  /// Communicator for the group
  MPI_Comm comm_;

  /// This process rank
  int process_rank_;

  /// First process in the group
  int process_first_;

  /// Last process (plus one) in the group 
  int process_last_plus_;

  /// Whether to use standard, buffered, synchronous, or ready sends
  send_type send_type_;

  /// Whether to use blocking sends
  bool send_blocking_;
  
  /// Whether to use blocking receives
  bool recv_blocking_;

};

#endif /* CONFIG_USE_MPI */

#endif /* PARALLEL_GROUP_PROCESS_MPI_HPP */

