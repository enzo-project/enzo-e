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
  /// @brief    Class for encapsulating different parallel technologies

public: // interface

  /// Initialize a Parallel object
  Parallel();

  /// Delete a Parallel object
  ~Parallel();

};

#endif /* PARALLEL_PARALLEL_HPP */

