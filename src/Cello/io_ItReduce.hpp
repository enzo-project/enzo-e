// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItReduce.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-26
/// @brief    [\ref Reduce] Base iterator class for reduction operators ItReduceMax, ItReduceSum, etc.

#ifndef IO_IT_REDUCE_HPP
#define IO_IT_REDUCE_HPP

class ItReduce  {

  /// @class    ItReduce
  /// @ingroup  Io
  /// @brief    [\ref Io] Iterator over Reduces in a Reduce

  
protected: //

  /// Prohibit creating new ItReduces directly: require ItReduce::create()
  ItReduce () throw ()
  { ERROR("ItReduce::ItReduce","ItReduce should not be called"); }

  /// Prohibit creating new ItReduces directly: require ItReduce::create()
  ItReduce (double value_start) throw ()
    : count_(0), value_(value_start) {}

public: // interface

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | count_;
    p | value_;
  }
#endif

  /// Static factory for creating ItReduce objects
  static ItReduce * create (reduce_enum reduce);

  /// Delete the ItReduce object
  virtual ~ItReduce () throw ()
  {} ;
  
  /// Reduce another value
  virtual void next (double value) throw() = 0;

  /// Reset the Iterator to the beginning
  virtual void first() throw() = 0;

  /// Return the current value of the reduction operator
  virtual double value() const throw() = 0;

protected: // attributes

  /// Number of values accumulated so far
  size_t count_;

  /// Current reduced value
  double value_;

};

#endif /* IO_IT_REDUCE_HPP */
