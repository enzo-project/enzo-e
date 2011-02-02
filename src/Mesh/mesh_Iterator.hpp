// $Id: mesh_Iterator.hpp 1942 2011-01-20 00:53:45Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_ITERATOR_HPP
#define MESH_ITERATOR_HPP

/// @file     mesh_Iterator.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 16:46:01 PST 2011
/// @brief    [\ref Mesh] Declaration of the Iterator abstract base class

class Iterator {

  /// @class    Iterator
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Abstract base class for all iterators

public: // interface

  /// Create an Iterator object
  Iterator () throw ()  {};

  /// Delete the Iterator object
  ~Iterator () throw () {};
  
  /// Iterate through all entities
  virtual void * operator++ () = 0;

};

#endif /* MESH_ITERATOR_HPP */
