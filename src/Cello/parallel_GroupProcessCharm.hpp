// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_GroupProcessCharm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 24 14:32:48 PDT 2011
/// @brief    [\ref Parallel] Interface for the GroupProcessCharm class

#ifndef PARALLEL_GROUP_PROCESS_CHARM_HPP
#define PARALLEL_GROUP_PROCESS_CHARM_HPP

#ifdef CONFIG_USE_CHARM

class GroupProcessCharm : public GroupProcess {

  /// @class    GroupProcessCharm
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] CHARM helper functions

public: // interface

  /// Initialize an GroupProcessCharm object
  GroupProcessCharm(int process_first     = 0,
		    int process_last_plus = -1) throw();

  /// Destructor
  virtual ~GroupProcessCharm() throw()
  {}

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    GroupProcess::pup(p);
    p | process_first_;
    p | process_last_plus_;
    // NOTE: change this function whenever attributes change
  }
#endif

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

private: // functions


private: // attributes

  /// First process in the group
  int process_first_;

  /// Last process (plus one) in the group 
  int process_last_plus_;

};

#endif /* CONFIG_USE_CHARM */

#endif /* PARALLEL_GROUP_PROCESS_CHARM_HPP */

