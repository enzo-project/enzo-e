// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItReduceSum.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-26
/// @brief    [\ref Reduce] Iterator class for averaging

#ifndef IO_IT_REDUCE_SUM_HPP
#define IO_IT_REDUCE_SUM_HPP

class ItReduceSum : public ItReduce {

  /// @class    ItReduceSum
  /// @ingroup  Io
  /// @brief    [\ref Io] Iterator over Reduces in a Reduce

public: //

  /// Prohibit creating new ItReduceSums directly: require ItReduceSum::create()
  ItReduceSum () throw ()
  : ItReduce(0.0)
  {}

public: // interface

  /// Delete the ItReduceSum object
  virtual ~ItReduceSum () throw ()
  {}
  
  /// Reduce another value
  virtual void next (double value) throw()
  { value_ += value; }

  /// Reset the Iterator to the beginning
  virtual void first() throw()
  { value_ = 0.0; }

  /// Return the current value of the reduction operator
  virtual double value() const throw()
  { return value_; }

private: // attributes

};

#endif /* IO_IT_REDUCE_SUM_HPP */
