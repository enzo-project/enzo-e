// $Id: field_FieldDescr.hpp 1261 2010-03-03 00:14:11Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELD_DESCR_HPP
#define FIELD_FIELD_DESCR_HPP

/// @file     field_FieldDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    Declaration for the FieldDescr class

enum field_action {
  field_action_unknown,  // Uninitialized action
  field_action_none,     // Do nothing if range exceeded
  field_action_set,      // Set field values to min / max if range exceeded
  field_action_warning,  // Issue warning if range exceeded
  field_action_error,    // Issue error if range exceeded
  field_action_timestep, // Retry with reduced timestep if range exceeded
  field_action_method    // Retry with alternate method if range exceeded
};

enum precision_type {
  precision_unknown,  //  unknown precision
  precision_default,  //  default precision, based on CONFIG_PRECISION_[SINGLE|DOUBLE]
  precision_32bit,    //  32-bit field data
  precision_64bit     //  64-bit field data
};

#include <string>
#include <memory>

class FieldDescr {

  /// @class    FieldDescr
  /// @ingroup  Field
  /// @brief    Interface for the FieldDescr class

public: // public

  /// Initialize a FieldDescr object
  FieldDescr(std::string name,
	     int dim = 3) throw();

  /// Set dimension
  void set_dimension (int dim)
  { if (1 <= dim && dim <= 3) {
      dim_ = dim;
    } else {
      ERROR_MESSAGE("FieldDescr::set_dimension","dim out of range");
    }
  };
      
  /// Return dimensionality
  int dimension () const throw() { return dim_;  };
      
  /// Return Field's name
  std::string name () const throw()
  { return name_; }

  /// Return centering of Field.  Assumes vector length is at least dim_
  const bool * centering () const throw()
  {
    return centering_;
  };

  /// Return centering of Field for the given axis
  void set_centering (int axis, bool value) throw()
  {
    if (0 <= axis && axis < dim_) 
      centering_[axis] = value; 
  };

  /// Return Field minimum value
  double min_value () const throw()
  { return min_value_; };

  /// Set Field minimum value
  void set_min_value (double value) throw()
  { min_value_ = value; };

  /// Return Field maximum value
  double max_value () const throw()
  { return max_value_; };

  /// Set Field maximum value
  void set_max_value (double value) throw()
  { max_value_ = value; };

  /// Return action on violating Field minimum value
  field_action min_action () const throw()
  { return min_action_; };

  /// Set action on violating Field minimum value
  void set_min_action (field_action action) throw()
  { min_action_ = action;  };

  /// Return action on violating Field maximum value
  field_action max_action () const throw()
  { return max_action_; };

  /// Set action on violating Field maximum value
  void set_max_action (field_action action) throw()
  { max_action_ = action;  };

  /// Return precision of Field
  precision_type precision () const throw()
  { return precision_; };

  /// Set precision of Field
  void set_precision (precision_type precision) throw()
  { 
    precision_ = (precision == precision_default) ? 
      default_precision_() : precision;
  };

private: // functions

  precision_type default_precision_ () {
#ifdef CONFIG_PRECISION_SINGLE
    return precision_32bit;
#endif
#ifdef CONFIG_PRECISION_DOUBLE
    return precision_64bit;
#endif
    return precision_unknown;
  };

private: // attributes

  /// Dimension of Field
  int dim_;

  /// String defining the field's name
  std::string name_;

  /// Cell centering
  bool * centering_;

  /// Minimum allowed value for the Field
  double min_value_;

  /// Maximum allowed value for the Field
  double max_value_;

  /// Action to perform if Field values go below min_
  field_action min_action_ ;

  /// Action to perform if Field values go above max_
  field_action max_action_ ;

  /// Precision of the Field data
  precision_type precision_;

};

#endif /* FIELD_FIELD_DESCR_HPP */

