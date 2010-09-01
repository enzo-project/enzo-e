// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef GLOBAL_GLOBAL_HPP
#define GLOBAL_GLOBAL_HPP

/// @file     global_Global.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file global_Global.hpp

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

#endif /* GLOBAL_GLOBAL_HPP */

