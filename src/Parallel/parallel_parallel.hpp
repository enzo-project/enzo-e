// $Id: parallel.hpp 1275 2010-03-09 01:05:14Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_PARALLEL_HPP
#define PARALLEL_PARALLEL_HPP

/// @file     parallel.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-16
/// @brief    Interface for the Parallel class

class Parallel {

  /// @class    Parallel
  /// @ingroup  Parallel
  /// @todo     Use singleton design pattern with variations: Mpi, Omp, etc.
  /// @brief    Class for encapsulating different parallel technologies

public: // interface

  /// Initialize
  virtual void initialize(int * argc = 0, char ***argv = 0) 
  {}

  /// Finalize
  virtual void finalize()
  {}

  /// Get size
  virtual int size()
  { return 1; }

  /// Get rank
  virtual int rank()
  { return 0; }

  virtual bool is_root()
  { return true; }

public: // static functions

  /// Get single instance of the Parallel object
  static Parallel * instance() throw ();

protected: // functions

  /// Initialize a Parallel object (singleton design pattern)
  Parallel() 
  {};

  /// Delete a Parallel object (singleton design pattern)
  ~Parallel() 
  {};

private: // attributes

  /// Single instance of the Parallel object (singleton design pattern)
  static Parallel * instance_;


};

#endif /* PARALLEL_PARALLEL_HPP */

