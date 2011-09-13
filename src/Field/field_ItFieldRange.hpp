// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_IT_FIELD_RANGE_HPP
#define FIELD_IT_FIELD_RANGE_HPP

/// @file     field_ItFieldRange.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief    [\ref Field] Declaration of the ItFieldRange iterator class

class ItFieldRange : public ItField {

  /// @class    ItFieldRange
  /// @ingroup  Field
  /// @brief    [\ref Field] Iterator over a range of Fields in a Block

public: // interface

  /// Create an iterator over integers 0 to count-1
  ItFieldRange ( size_t first, size_t last ) throw ()
    : ItField(0), first_(first), last_(last)
  { }

  /// Go to the first value
  virtual void first () throw()
  { index1_ = first_; }

  /// Go to the next value
  virtual void next () throw()
  { if (index1_ <= last_) index1_++; }

  /// Return the current value
  virtual int value () const throw()
  { return int(index1_); }

  /// Return whether iterating is complete
  virtual bool done () const throw()
  { return (index1_ > last_); }

private: // attributes

  /// First value
  size_t first_;

  /// Last value
  size_t last_;

};

#endif /* FIELD_IT_FIELD_RANGE_HPP */
