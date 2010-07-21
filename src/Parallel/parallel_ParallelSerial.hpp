// $Id: parallel.hpp 1275 2010-03-09 01:05:14Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_ParallelSerial.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Rename Parallel  ParallelSerial
/// @todo     Add ParallelCharm (for CkExit(), etc.)
/// @date     2009-10-16
/// @brief    Interface for the Parallel class

#ifndef PARALLEL_PARALLEL_SERIAL_HPP
#define PARALLEL_PARALLEL_SERIAL_HPP

class ParallelSerial : public Parallel {

  /// @class    ParallelSerial
  /// @ingroup  Parallel
  /// @todo     Split into ParallelProcesses and ParallelThreads or similar
  /// @brief    Class for encapsulating different, possibly multiple,
  /// parallel technologies

public: // interface

  /// Initialize a ParallelSerial object (singleton design pattern)
  ParallelSerial() 
    : Parallel()
  {};

  /// Initialize
  virtual void initialize(int * argc = 0, char ***argv = 0) 
  { 
      printf ("ParallelSerial::initialize() this = %p\n",this);
      set_initialized_(true);
  }

  /// Finalize
  virtual void finalize()
  { set_initialized_(false); }

  /// Abort execution abruptly
  virtual void abort()
  { exit (1); }

  /// Exit the program
  virtual void halt()
  { exit (0); }

  /// Get total number of processors
  virtual int process_count()
  { return 1; }

  /// Get rank of this process
  virtual int process_rank()
  { return 0; }

  /// Get total number of threads in this node
  virtual int thread_count()
  { return 1; }

  /// Get rank of this thread
  virtual int thread_rank()
  { return 0; }

  /// Get rank
  virtual std::string name()
  { return "0"; }

  /// Return whether this is the root process
  virtual bool is_root()
  { return true; }

};

#endif /* PARALLEL_PARALLEL_SERIAL_HPP */

