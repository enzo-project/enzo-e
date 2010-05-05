// $Id: strict_auto_ptr.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef STRICT_AUTO_PTR_HPP
#define STRICT_AUTO_PTR_HPP

/// @file     strict_auto_ptr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed May  5 13:50:32 PDT 2010
/// @brief    Implementation of a strict auto_ptr.  See C++ FAQs 2nd ed p.427

#include <memory>

template<class T>
class strict_auto_ptr : public std::auto_ptr<T> {
 public:
  strict_auto_ptr(T* p = NULL) throw() : std::auto_ptr<T>(p) { }
 private:
  strict_auto_ptr (const strict_auto_ptr&) throw();
  void operator = ( const strict_auto_ptr&) throw();
};

#endif /* STRICT_AUTO_PTR_HPP */


