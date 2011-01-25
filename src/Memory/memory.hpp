// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MEMORY_HPP
#define MEMORY_HPP

/// @file     memory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 11 17:20:03 PST 2010
/// @brief    Include file for the \ref Memory component 

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stack>
#include <memory>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "error.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "memory_Memory.hpp"

/// strict_auto_ptr class template
template<class T>
class strict_auto_ptr : public std::auto_ptr<T> {
public:
  strict_auto_ptr(T* p = NULL) throw() : std::auto_ptr<T>(p) { }
private:
  strict_auto_ptr (const strict_auto_ptr&) throw();
  void operator = ( const strict_auto_ptr&) throw();
};

#endif /* MEMORY_HPP */

