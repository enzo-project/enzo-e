// $Id: parallel.hpp 1275 2010-03-09 01:05:14Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Parallel.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Add ParallelCharm (for CkExit(), etc.)
/// @date     2009-10-16
/// @brief    Interface for the Parallel class

#ifndef PARALLEL_PARALLEL_HPP
#define PARALLEL_PARALLEL_HPP

#include <string>
#include <stdio.h>
#include <stdlib.h>

class Parallel {

  /// @class    Parallel
  /// @ingroup  Parallel
  /// @brief    Base class for Parallel objects, e.g. ParallelMpi, ParallelSerial

public: // interface

  /// Initialize a Parallel object (singleton design pattern)
  Parallel() :
    initialized_(false)
  {};

  /// Initialize
  virtual void initialize(int * argc = 0, char ***argv = 0) = 0;

  /// Finalize
  virtual void finalize() = 0;

  /// Abort execution abruptly
  virtual void abort() = 0;

  /// Exit the program
  virtual void halt() = 0;

  /// Get total number of processors
  virtual int process_count() = 0;

  /// Get rank of this process
  virtual int process_rank() = 0;

  /// Get total number of threads in this node
  virtual int thread_count() = 0;

  /// Get rank of this thread
  virtual int thread_rank() = 0;

  /// Get rank
  virtual std::string name() = 0;

  /// Return whether this is the root process
  virtual bool is_root() = 0;

protected:

  /// Set whether class is initialized
  void set_initialized_(bool initialized)
  {initialized_ = initialized; };

  /// Return whether class is initialized
  bool is_initialized_()
  { return initialized_ ; };

private: // attributes

  bool initialized_;

};

#endif /* PARALLEL_PARALLEL_HPP */

