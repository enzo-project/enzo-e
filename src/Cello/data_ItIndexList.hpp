// See LICENSE_CELLO file for license and copyright information

/// @file     data_ItIndexList.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief    [\ref Data] Declaration of the ItIndexList iterator class

#ifndef DATA_IT_INDEX_LIST_HPP
#define DATA_IT_INDEX_LIST_HPP

class ItIndexList : public ItIndex {

  /// @class    ItIndexList
  /// @ingroup  Data
  /// @brief    [\ref Data] Iterator over a list of indices

public: // interface

  /// Create an iterator over integers 0 to count-1
  ItIndexList () throw ()
  : ItIndex (),
    index_(0),
    values_()
  { }

  /// Virtual destructor
  virtual ~ItIndexList () throw()
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(ItIndexList);

  /// Charm++ PUP::able migration constructor
  ItIndexList (CkMigrateMessage *m)
    : ItIndex (m),
      index_(0),
      values_()
  { }
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    // NOTE: change this function whenever attributes change
    ItIndex::pup(p);
    p | index_;
    p | values_;
  }

  /// Append a value to the list of values
  void append (int value) 
  { values_.push_back(value); 
    size_++;}

  /// Go to the first value
  virtual void first () throw()
  { index_ = 0; }

  /// Go to the next value
  virtual void next () throw()
  { if (index_ < values_.size()) index_++; }

  /// Return the current value.  Should not be called if done() == true
  virtual int value () const throw()
  { return values_[index_]; }

  /// Return whether iterating is complete
  virtual bool done () const throw()
  { return (index_ >= values_.size()); }

private: // attributes

  /// Index of the current value
  size_t index_;

  /// List of values
  std::vector <int> values_;
};

#endif /* DATA_IT_INDEX_LIST_HPP */
