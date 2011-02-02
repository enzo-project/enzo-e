// $Id: method_Iterator.hpp 1942 2011-01-20 00:53:45Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef METHOD_ITERATOR_HPP
#define METHOD_ITERATOR_HPP

/// @file     method_Iterator.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 16:46:01 PST 2011
/// @brief    [\ref Method] Declaration of the Iterator abstract base class

class Iterator {

  /// @class    Iterator
  /// @ingroup  Method
  /// @brief    [\ref Method] Abstract base class for all iterators

public: // interface

  /// Create an Iterator object
  Iterator () throw ()  {};

  /// Delete the Iterator object
  ~Iterator () throw () {};
  
  /// Iterate through all entities
  virtual void * operator++ () = 0;

};

#endif /* METHOD_ITERATOR_HPP */
