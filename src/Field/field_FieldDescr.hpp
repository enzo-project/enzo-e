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
  FieldDescr(int dim) throw();

  /// Add a new field to the list of known fields
  int add_field (std::string name) throw()
  { 
    int index = field_count();
    name_.push_back(name); 
    bool * centering = new bool [dim_];
    for (int i=0; i < dim_; i++) centering[i] = true;
    centering_. push_back(centering);
    min_value_. push_back(0);
    max_value_. push_back(0);
    min_action_.push_back(field_action_none);
    max_action_.push_back(field_action_none);
    precision_. push_back(default_precision_());
    return index;
  };

  /// Return the number of Fields
  int field_count() const throw()
  { return name_.size(); }

  /// Return indexed Field's name
  std::string name (int i) const throw()
  { return (0 <= i && i < field_count()) ? 
      name_[i] : ""; };

  /// Return named Field's index
  int index (std::string name) const throw()
  { 
    for (int index=0; index < field_count(); index++) {
      if (this->name(index) == name) return index;
    }
    return -1; // Uh oh
  }

  /// Return centering of Field.  Assumes vector length is at least dim_
  const bool * centering (int index) const throw()
  {
    return (0 <= index && index < field_count()) ? 
      centering_[index] : 0; 
  };

  /// Return centering of Field for the given axis.  Assumes index and axis are in range.
  void set_centering (int index,int axis, bool value) throw()
  {
    if (0 <= index && index < field_count() && 0 <= axis && axis < dim_) 
      centering_[index][axis] = value; 
  };

  /// Return ith Field minimum value
  double min_value (int index) const throw()
  { return (0 <= index && index < field_count()) ? 
      min_value_[index] : 0.0; };

  /// Set ith Field minimum value
  void set_min_value (int index,double value) throw()
  { 
    if (0 <= index && index < field_count()) {
      min_value_[index] = value;
    }
  };

  /// Return ith Field maximum value
  double max_value (int index) const throw()
  { return (0 <= index && index < field_count()) ? 
      max_value_[index] : 0.0; };

  /// Set ith Field maximum value
  void set_max_value (int index,double value) throw()
  { 
    if (0 <= index && index < field_count()) {
      max_value_[index] = value;
    }
  };

  /// Return action on violating ith Field minimum value
  field_action min_action (int index) const throw()
  { return (0 <= index && index < field_count()) ? 
      min_action_[index] : field_action_unknown; };

  /// Set action on violating ith Field minimum value
  void set_min_action (int index,field_action action) throw()
  { 
    if (0 <= index && index < field_count()) {
      min_action_[index] = action;
    }
  };

  /// Return action on violating ith Field maximum value
  field_action max_action (int index) const throw()
  { return (0 <= index && index < field_count()) ? 
      max_action_[index] : field_action_unknown; };

  /// Set action on violating ith Field maximum value
  void set_max_action (int index,field_action action) throw()
  { 
    if (0 <= index && index < field_count()) {
      max_action_[index] = action;
    }
  };

  /// Return precision of ith Field
  precision_type precision (int index) const throw()
  {
    return (0 <= index && index < field_count()) ? 
      precision_[index] : precision_unknown; };

  /// Set precision of ith Field
  void set_precision (int index, precision_type precision) throw()
  { 
    if (precision == precision_default) {
      precision = default_precision_();
    }
    if (0 <= index && index < field_count()) {
      precision_[index] = precision;
    }
  };

  /// 
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

  /// Dimension of Fields
  int dim_;

  /// String defining the field's name
  std::vector < std::string> name_;

  /// Cell centering
  std::vector < bool *> centering_;

  /// Minimum allowed value for the Field
  std::vector < double> min_value_;

  /// Maximum allowed value for the Field
  std::vector < double> max_value_;

  /// Action to perform if Field values go below min_
  std::vector < field_action> min_action_ ;

  /// Action to perform if Field values go above max_
  std::vector < field_action> max_action_ ;

  /// Precision of the Field data
  std::vector < precision_type> precision_;

};

#endif /* FIELD_FIELD_DESCR_HPP */

