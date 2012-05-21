// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_CountersMpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-19
/// @brief    [\ref Lcaperf] Declaration of the CountersMpi class

#ifndef LCAPERF_COUNTERS_MPI_HPP
#define LCAPERF_COUNTERS_MPI_HPP

#ifdef CONFIG_USE_MPI

namespace lca {

class CountersMpi : public CountersUser {

  /// @class    CountersMpi
  /// @ingroup  lcaperf
  /// @brief    [\ref Lcaperf] Mpi counters

public: // interface

  /// Constructor
  CountersMpi() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~CountersMpi() throw();

  /// Copy constructor
  CountersMpi(const CountersMpi & counters) throw();

  /// Assignment operator
  CountersMpi & operator= (const CountersMpi & counters) throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    //    p |  ip_mpi_;
    //    p |  np_mpi_;
    //    p |  ip_node_;
    //    p |  np_node_;
    ERROR("CountersMpi::pup",
	  "Charm++ pup() function called in MPI configuration");

  }
#endif

  //----------------------------------------------------------------------

  /// Initialize CountersMpi 
  virtual void initialize();

  /// Update send-related counters
  inline static void count_send(int n, MPI_Datatype t, int c=1)
  {
    int s;
    MPI_Type_size(t,&s);

    mpi_send_time_  += mpi_call_time_;
    mpi_send_bytes_ += c*n*s;
    mpi_send_calls_ += c;
  }

  //----------------------------------------------------------------------

  /// Update recv-related counters
  inline static void count_recv(int n, MPI_Datatype t, int c=1)
  {
    int s;
    MPI_Type_size(t,&s);

    mpi_recv_time_  += mpi_call_time_;
    mpi_recv_bytes_ += c*n*s;
    mpi_recv_calls_ += c;
  }

  //----------------------------------------------------------------------

  /// Update sync-related counters
  inline static void count_sync(MPI_Comm c)
  {
    int np;
    MPI_Comm_size(c,&np);

    mpi_sync_time_  += mpi_call_time_;
    mpi_sync_procs_ += np;
    mpi_sync_calls_ += 1;
  }

  //----------------------------------------------------------------------
  /// Start timing the current MPI function
  inline static void timer_start ()
  { 
    mpi_time_start_ = MPI_Wtime();
  };

  //----------------------------------------------------------------------
  /// Start timing the current MPI function
  inline static void timer_stop ()
  {
    // Convert time from s to us
    mpi_call_time_ = 1000000L * (MPI_Wtime() - mpi_time_start_);
    mpi_time_  += mpi_call_time_;
    mpi_calls_ += 1;
  };



protected: // functions

  /// Create and start a new set of counters for the current key
  virtual long long * start_();

  /// stop a set of counters for the current key
  virtual void stop_(long long * );

  /// Function for copying static MPI counters to lcaperf user counters
  void assign_all_();

  //----------------------------------------------------------------------

private: // functions

  //----------------------------------------------------------------------

private: // attributes

  // Rank of this MPI process
  int ip_mpi_;

  // Number of MPI processes
  int np_mpi_;

  // rank of this "node" out of all nodes
  int ip_node_;

  // Number of nodes
  int np_node_;

private: // static attributes

  /// Number of calls to MPI functions
  static long long mpi_calls_;

  /// Time spent in all MPI calls
  static long long mpi_time_;

  /// Time spent in sends
  static long long mpi_send_time_;

  /// Number of bytes sent
  static long long mpi_send_bytes_;

  /// Number of message packets sent
  static long long mpi_send_calls_;

  /// Time spent in receives
  static long long mpi_recv_time_;

  /// Number of bytes received
  static long long mpi_recv_bytes_;

  /// Number of message packets received
  static long long mpi_recv_calls_;

  /// Time spent in barrier calls
  static long long mpi_sync_time_;

  /// Number of processes involved in all barrier calls
  static long long mpi_sync_procs_;

  /// Number of barrier calls
  static long long mpi_sync_calls_;

  /// Start time in current MPI call: used to compute mpi_time_
  static double mpi_time_start_;

  /// Time spent in latest MPI call
  static long long mpi_call_time_;

};

}

#endif /* CONFIG_USE_MPI */

#endif /* LCAPERF_COUNTERS_MPI_HPP */

