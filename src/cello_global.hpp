// $Id: cello_global.hpp 1746 2010-09-01 00:29:25Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef CELLO_GLOBAL_HPP
#define CELLO_GLOBAL_HPP

/// @file     cello_global.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Container for global objects

#include "error.hpp"
#include "memory.hpp"
#include "monitor.hpp"
#include "parameters.hpp"
#include "parallel.hpp"

class Global {

  /// @class    Global
  /// @ingroup  Global
  /// @brief    Object for storing cross-cutting objects, including monitor, error, parameters

public: // interface

  /// Constructor
  Global() throw()
  {
    error_      = new Error;
    monitor_    = new Monitor();
    parameters_ = new Parameters(monitor_);
    //    memory_     = new Memory;
  }

  /// Destructor
  ~Global() throw()
  {
    delete error_;
    delete parameters_;
    delete monitor_;
    //    delete [] memory_;
  }

private: // prohibit copy constructor

  /// Copy constructor
  Global(const Global & global) throw()
  {
  }

private: // prohibit assignment

  /// Assignment operator
  Global & operator= (const Global & global) throw();

public:

  /// Access error object

  Error * error() { return error_; };

  /// Access monitor object

  Monitor * monitor() { return monitor_; };

  /// Access parameters object

  Parameters * parameters() { return parameters_; };

//   /// Access memory object

//   Memory * memory() { return memory_; };


private: // attributes

  /// Error object
  Error * error_;

  /// Monitor object
  Monitor * monitor_;

  /// Parameters object
  Parameters * parameters_;

//   /// Memory object
//   Memory * memory_;

  ///  Performance object
  // Performance * performance_;

};

#endif /* CELLO_GLOBAL_HPP */

