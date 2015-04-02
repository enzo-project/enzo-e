// See LICENSE_CELLO file for license and copyright information

/// @file     data_ItField.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief [\ref Data] Declaration of the ItField abstract iterator
/// base class

#ifndef DATA_IT_FIELD_HPP
#define DATA_IT_FIELD_HPP

class ItField : public PUP::able 
{

  /// @class    ItField
  /// @ingroup  Field
  /// @brief    [\ref Data] Abstract iterator base class for Field indices

public: // interface

  /// Create an iterator over integers 0 to count-1
  ItField () throw ()
  { }

  /// Virtual destructor
  virtual ~ItField () throw()
  { }

  /// Charm++ PUP::able declarations
  PUPable_abstract(ItField);

  /// Charm++ PUP::able migration constructor
  ItField (CkMigrateMessage *m) : PUP::able(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  {
    TRACEPUP;
    PUP::able::pup(p); 
  };

  /// Go to the first value
  virtual void first () throw() = 0;

  /// Go to the next value
  virtual void next () throw() = 0;

  /// Return the current value
  virtual int value () const throw() = 0;

  /// Return whether iterating is complete
  virtual bool done () const throw() = 0;

protected: // attributes

};

#endif /* DATA_IT_FIELD_HPP */
