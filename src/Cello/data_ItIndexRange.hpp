// See LICENSE_CELLO file for license and copyright information

/// @file     data_ItIndexRange.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief    [\ref Data] Declaration of the ItIndexRange iterator class

#ifndef DATA_IT_INDEX_RANGE_HPP
#define DATA_IT_INDEX_RANGE_HPP

class ItIndexRange : public ItIndex {

  /// @class    ItIndexRange
  /// @ingroup  Data
  /// @brief    [\ref Data] Iterator over a range of Indices

public: // interface

  /// Empty constructor for Charm++ pup()
  ItIndexRange() throw() : ItIndex() {}
  
  /// Create an iterator over integers first to last
  ItIndexRange ( size_t first, size_t last ) throw ()
    : ItIndex(), first_(first), last_(last)
  { size_ = last_ - first_ + 1; }

  /// Create an iterator over integers 0 to count - 1first to last
  ItIndexRange ( size_t count ) throw ()
    : ItIndex(), index_(0), first_(0), last_(count - 1)
  { size_ = 0;}

  /// Virtual destructor
  virtual ~ItIndexRange () throw ()
  {}

  /// Charm++ PUP::able declarations
  PUPable_decl(ItIndexRange);

  /// Charm++ PUP::able migration constructor
  ItIndexRange (CkMigrateMessage *m)
    : ItIndex (m),
      index_(0),
      first_(0),
      last_(0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    // NOTE: change this function whenever attributes change
    ItIndex::pup(p);
    p | index_;
    p | first_;
    p | last_;
  }

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

#endif /* DATA_IT_INDEX_RANGE_HPP */
