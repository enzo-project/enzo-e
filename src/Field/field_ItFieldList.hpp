// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_IT_FIELD_LIST_HPP
#define FIELD_IT_FIELD_LIST_HPP

/// @file     field_ItFieldList.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief    [\ref Field] Declaration of the ItFieldList iterator class

class ItFieldList : public ItField {

  /// @class    ItFieldList
  /// @ingroup  Field
  /// @brief    [\ref Field] Iterator over Fields in a Block

public: // interface

  /// Create an iterator over integers 0 to count-1
  ItFieldList ( size_t count ) throw ()
    : ItField (count), values_()
  { values_.reserve(count); }

  /// Append a value to the list of values
  void append (int value) 
  { values_.push_back(value); }

  /// Go to the next value
  virtual void next () throw()
  { if (index1_ < values_.size()) index1_++; }

  /// Return the current value
  virtual int value () const throw()
  { return values_[index1_]; }

  /// Return whether iterating is complete
  virtual bool done () const throw()
  { return (index1_ >= values_.size()); }

  private: // attributes

  std::vector <int> values_;
};

#endif /* FIELD_IT_FIELD_LIST_HPP */
