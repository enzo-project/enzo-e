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

public: //

  /// Prohibit creating new ItReduceAvgs directly: require ItReduceAvg::create()
  ItReduceAvg () throw ();

public: // interface

  /// Delete the ItReduceAvg object
  virtual ~ItReduceAvg () throw ();
  
  /// Reduce another value
  virtual void next (double value) throw();

  /// Reset the Iterator to the beginning
  virtual void first() throw();

  /// Return the current value of the reduction operator
  virtual double value() const throw();

private: // attributes

};

#endif /* IO_IT_REDUCE_AVG_HPP */
