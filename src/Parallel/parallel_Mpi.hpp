// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_mpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Oct 15 10:40:37 PDT 2009 
/// @brief    Interface for the ParallelMpi class

#ifndef PARALLEL_MPI_HPP
#define PARALLEL_MPI_HPP

enum enum_send_type {
  send_type_standard,
  send_type_buffered,
  send_type_synchronous,
  send_type_ready };

enum enum_data_type {
  data_type_single,
  data_type_double };

class ParallelMpi : public Parallel {

  /// @class    ParallelMpi
  /// @ingroup  Parallel
  /// @brief    MPI helper functions

public: // interface

  /// Initialize MPI (virtual)
  void initialize(int * argc, char ***argv);

  /// Finalize MPI
  void finalize();

  /// Abort execution abruptly
  void abort();

  /// Exit the program
  void halt();

  /// Get MPI size
  int process_count() 
  { return size_; }

  /// Get MPI rank
  int process_rank() 
  { return rank_; }

  virtual std::string name()
  { return name_; }

public: // interface

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

public: // static functions

  /// Get single instance of the Parallel object
  static ParallelMpi * instance() throw ()
  { return & ParallelMpi::instance_mpi_; }

public: // virtual

  /// Return whether this is the root process
  virtual bool is_root()
  { return rank_ == 0; }

protected: // functions

  /// Initialize an ParallelMpi object (singleton design pattern)
  ParallelMpi()
    : Parallel(),
      size_(1),
      rank_(0),
      name_("0"),
      send_blocking_(true),
      recv_blocking_(true),
      send_type_(send_type_standard)
  {};

private: // static functions

  /// Single instance of the Parallel object (singleton design pattern)
  static ParallelMpi instance_mpi_;

private: // attributes

  /// MPI global size
  int size_;

  /// MPI global rank
  int rank_;

  /// MPI global name
  std::string name_;

  /// Whether to use blocking sends
  bool send_blocking_;
  
  /// Whether to use blocking receives
  bool recv_blocking_;

  /// Whether to use standard, buffered, synchronous, or ready sends
  enum_send_type send_type_;

protected: // static attributes

};

#endif /* PARALLEL_MPI_HPP */

