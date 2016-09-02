// See LICENSE_CELLO file for license and copyright information

/// @file     data_ItIndex.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief [\ref Data] Declaration of the ItIndex abstract iterator
/// base class

#ifndef DATA_IT_INDEX_HPP
#define DATA_IT_INDEX_HPP

class ItIndex : public PUP::able 
{

  /// @class    ItIndex
  /// @ingroup  Data
  /// @brief    [\ref Data] Abstract iterator base class for Index indices

public: // interface

  /// Create an iterator over integers 0 to count-1
  ItIndex () throw () : size_(0)
  { }

  /// Virtual destructor
  virtual ~ItIndex () throw()
  { }

  /// Charm++ PUP::able declarations
  PUPable_abstract(ItIndex);

  /// Charm++ PUP::able migration constructor
  ItIndex (CkMigrateMessage *m)
    : PUP::able(m),
      size_(0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  {
    TRACEPUP;
    PUP::able::pup(p); 
    p | size_;
  };

  /// Number of elements contained in the list
  int size () const throw() { return size_; }

  /// Go to the first value
  virtual void first () throw() = 0;

  /// Go to the next value
  virtual void next () throw() = 0;

  /// Return the current value
  virtual int value () const throw() = 0;

  /// Return whether iterating is complete
  virtual bool done () const throw() = 0;

protected: // attributes

  int size_;
};

#endif /* DATA_IT_INDEX_HPP */
