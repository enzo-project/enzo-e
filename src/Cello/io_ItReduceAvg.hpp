// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItReduceAvg.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-26
/// @brief    [\ref Reduce] Iterator class for averaging

#ifndef IO_IT_REDUCE_AVG_HPP
#define IO_IT_REDUCE_AVG_HPP

class ItReduceAvg : public ItReduce {

  /// @class    ItReduceAvg
  /// @ingroup  Io
  /// @brief    [\ref Io] Iterator over Reduces in a Reduce

public: // interface

  /// Create an ItReduceAvg object
  ItReduceAvg () throw ()
  : ItReduce(0.0)
  { }

  /// Delete the ItReduceAvg object
  virtual ~ItReduceAvg () throw ()
  { }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    ItReduce::pup(p);
  }
  
  /// Reduce another value
  virtual void next (double value) throw()
  {
    count_++;
    value_ += value;
  }

  /// Reset the Iterator to the beginning
  virtual void first() throw()
  {
    count_=0;
    value_ = 0.0;
  }

  /// Return the current value of the reduction operator
  virtual double value() const throw()
  { return value_ / count_; }

private: // attributes

};

#endif /* IO_IT_REDUCE_AVG_HPP */
