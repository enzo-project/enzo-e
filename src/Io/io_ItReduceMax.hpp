// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItReduceMax.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-26
/// @brief    [\ref Reduce] Iterator class for averaging

#ifndef IO_IT_REDUCE_MAX_HPP
#define IO_IT_REDUCE_MAX_HPP

class ItReduceMax : public ItReduce {

  /// @class    ItReduceMax
  /// @ingroup  Io
  /// @brief    [\ref Io] Iterator over Reduces in a Reduce

public: //

  /// Prohibit creating new ItReduceMaxs directly: require
  /// ItReduceMax::create()
  ItReduceMax () throw ()
  : ItReduce(std::numeric_limits<double>::min())
  { }


public: // interface

  /// Delete the ItReduceMax object
  virtual ~ItReduceMax () throw ()
  { }
  
  /// Reduce another value
  virtual void next (double value) throw()
  { value_ = MAX(value,value_); }

  /// Reset the Iterator to the beginning
  virtual void first() throw()
  { value_ = std::numeric_limits<double>::min(); }


  /// Return the current value of the reduction operator
  virtual double value() const throw()
  { return value_; }

private: // attributes

};

#endif /* IO_IT_REDUCE_MAX_HPP */
