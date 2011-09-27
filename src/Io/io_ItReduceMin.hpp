// See LICENSE_CELLO file for license and copyright information

#ifndef IO_IT_REDUCE_MIN_HPP
#define IO_IT_REDUCE_MIN_HPP

/// @file     io_ItReduceMin.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-26
/// @brief    [\ref Reduce] Iterator class for averaging


class ItReduceMin : public ItReduce {

  /// @class    ItReduceMin
  /// @ingroup  Io
  /// @brief    [\ref Io] Iterator over Reduces in a Reduce

public: //

  /// Prohibit creating new ItReduceMins directly: require ItReduceMin::create()
  ItReduceMin () throw ();

public: // interface

  /// Delete the ItReduceMin object
  virtual ~ItReduceMin () throw ();
  
  /// Reduce another value
  virtual void next (double value) throw();

  /// Reset the Iterator to the beginning
  virtual void first() throw();

  /// Return the current value of the reduction operator
  virtual double value() const throw();

private: // attributes

};

#endif /* IO_IT_REDUCE_MIN_HPP */
