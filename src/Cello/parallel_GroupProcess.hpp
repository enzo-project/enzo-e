// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_GroupProcess.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul 22 12:36:38 PDT 2010
/// @brief    [\ref Parallel] Declaration of the GroupProcess class

#ifndef PARALLEL_GROUP_PROCESS_HPP
#define PARALLEL_GROUP_PROCESS_HPP

class Reduce;

class GroupProcess
{

  /// @class    GroupProcess
  /// @ingroup  Parallel  
  /// @brief    [\ref Parallel] Group of distributed memory processes

 public: // static interface

  static GroupProcess * create (int process_first     = 0,
				int process_last_plus = -1) throw();

 protected: // interface

  /// Protected since GroupProcess objects must be created with create()
  GroupProcess(int size = 1, int rank = 0) throw()
    : size_(size),
    rank_(rank)
    {  }

 public: // interface

  /// Destructor
  virtual ~GroupProcess() throw()
  {}

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | size_;
    p | rank_;
  }
#endif

  //--------------------------------------------------
  // Shared with GroupThread, but removed from deleted common "Group"
  // base class because of name-clash with CHARM++
  //--------------------------------------------------

  /// Number of compute elements in the GroupProcess
  int size() const throw()
  { return size_; };

  /// Rank of the compute element in the GroupProcess
  int rank() const throw()
  {  return rank_; };

  /// True iff rank() is 0
  bool is_root() const throw()
  {  return rank_==0; };

  /// Synchronize between all compute elements in the GroupProcess
  virtual void barrier() const throw() { };

  /// Synchronize between two compute elements in the GroupProcess
  virtual void sync(int rank, int tag=0) const throw() { };

  //--------------------------------------------------

  /// Initiate sending an array
  virtual void * send_begin
    (int rank_dest, void * buffer, int size, int tag=0) const throw() = 0;

  /// Clean up after sending an array
  virtual void send_end(void * handle) const throw() = 0;

  /// Initiate receiving an array
  virtual void * recv_begin
    (int rank_source, void * buffer, int size, int tag=0) const throw() = 0;

  /// Clean up after receiving an array
  virtual void recv_end(void * handle) const throw() = 0;

  /// Send and receive data between two processes
  virtual void send_recv
  (int rank, void * buffer, int size, int tag=0) const throw() = 0;

  /// Wait until the send or receive operation completes
  virtual void wait(void * handle) const throw() = 0;

  /// Test completeness of sending or receiving an array
  virtual bool test (void * handle) const throw() = 0;

  //--------------------------------------------------

  /// Add an array to a list of arrays to send in bulk
  virtual void bulk_send_add
  (int rank_dest, void * buffer, int size, int tag=0) const throw() = 0;

  /// Initiate a bulk send of multiple arrays
  virtual void * bulk_send() const throw() = 0;

  /// Add an array to a list of arrays to receive in bulk
  virtual void bulk_recv_add
  (int rank_source, void * buffer, int size, int tag=0) const throw() = 0;

  /// Initiate a bulk receive of multiple arrays
  virtual void * bulk_recv() const throw() = 0;

  /// Complete a bulk send or receive of multiple arrays
  virtual void bulk_wait(void * handle) const throw() = 0;

  //--------------------------------------------------

  /// Create a Reduce object for this GroupProcess
  virtual Reduce * create_reduce () const throw ()= 0;

 protected: // attributes

  /// Number of compute elements in the GroupProcess
  int size_;

  /// Rank of this compute element in the GroupProcess
  int rank_;

};

#endif /* PARALLEL_GROUP_PROCESS_HPP */

