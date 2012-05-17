// See LICENSE_CELLO file for license and copyright information

/// @file     field_ItField.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief [\ref Field] Declaration of the ItField abstract iterator
/// base class

#ifndef FIELD_IT_FIELD_HPP
#define FIELD_IT_FIELD_HPP

class ItField {

  /// @class    ItField
  /// @ingroup  Field
  /// @brief    [\ref Field] Abstract iterator base class for Field indices

public: // interface

  /// Create an iterator over integers 0 to count-1
  ItField () throw ()
  { }

  /// Virtual destructor
  virtual ~ItField () throw()
  { }

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change

    // this space intentionally blank
  }
#endif

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

#endif /* FIELD_IT_FIELD_HPP */
