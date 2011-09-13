// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_IT_FIELD_HPP
#define FIELD_IT_FIELD_HPP

/// @file     field_ItField.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief    [\ref Field] Declaration of the ItField iterator class

class ItField {

  /// @class    ItField
  /// @ingroup  Field
  /// @brief    [\ref Field] Iterator over Fields in a Block

public: // interface

  /// Create an iterator over integers 0 to count-1
  ItField ( size_t count ) throw ()
    : index1_(0), count_ (count)
  { }

  /// Go to the first value
  virtual void first () throw()
  { index1_ = 0; }

  /// Go to the next value
  virtual void next () throw()
  { if (index1_ < count_) index1_++; }

  /// Return the current value
  virtual int value () const throw()
  { return int(index1_); }

  /// Return whether iterating is complete
  virtual bool done () const throw()
  { return (index1_ >= count_); }

protected: // attributes

  /// Current index plus one

  size_t index1_;

private: // attributes

  /// Maximum count
  size_t count_;

};

#endif /* FIELD_IT_FIELD_HPP */
