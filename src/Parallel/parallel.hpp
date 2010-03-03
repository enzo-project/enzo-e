// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

/// @file     parallel.hpp
/// @todo     Split into separate component include and class define files
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-16
/// @brief    Interface for the Parallel class and component include file

class Parallel {

  /// @class    Parallel
  /// @ingroup  Parallel
  /// @brief    Class for encapsulating different parallel technologies

public: // interface

  /// Initialize a Parallel object
  Parallel();

  /// Delete a Parallel object
  ~Parallel();

};

#include "mpi.hpp"

#endif /* PARALLEL_HPP */

