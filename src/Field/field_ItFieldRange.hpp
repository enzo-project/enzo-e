// See LICENSE_CELLO file for license and copyright information

/// @file     field_ItFieldRange.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief    [\ref Field] Declaration of the ItFieldRange iterator class

#ifndef FIELD_IT_FIELD_RANGE_HPP
#define FIELD_IT_FIELD_RANGE_HPP

class ItFieldRange : public ItField {

  /// @class    ItFieldRange
  /// @ingroup  Field
  /// @brief    [\ref Field] Iterator over a range of Fields in a Block

public: // interface

  /// Create an iterator over integers first to last
  ItFieldRange ( size_t first, size_t last ) throw ()
    : ItField(), first_(first), last_(last)
  { }

  /// Create an iterator over integers 0 to count - 1first to last
  ItFieldRange ( size_t count ) throw ()
    : ItField(), first_(0), last_(count - 1)
  { }

  /// Go to the first value
  virtual void first () throw()
  { index_ = first_; }

  /// Go to the next value
  virtual void next () throw()
  { if (index_ <= last_) index_++; }

  /// Return the current value
  virtual int value () const throw()
  { return int(index_); }

  /// Return whether iterating is complete
  virtual bool done () const throw()
  { return (index_ > last_); }

private: // attributes

  /// Index of the current value
  size_t index_;

  /// First value
  size_t first_;

  /// Last value
  size_t last_;

};

#endif /* FIELD_IT_FIELD_RANGE_HPP */
