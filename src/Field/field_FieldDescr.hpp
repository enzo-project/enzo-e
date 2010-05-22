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
  field_action_assign,   // Assign field values to min / max if range exceeded
  field_action_warning,  // Issue warning if range exceeded
  field_action_error,    // Issue error if range exceeded
  field_action_timestep, // Retry with reduced timestep if range exceeded
  field_action_method    // Retry with alternate method if range exceeded
};

enum precision_type {
  precision_unknown,   //  unknown precision
  precision_default,   //  default precision, based on CONFIG_PRECISION_[SINGLE|DOUBLE]
  precision_half,      //  16-bit field data
  precision_single,    //  32-bit field data
  precision_double,    //  64-bit field data
  precision_extended,  //  80-bit (ala Intel) field data
  pricision_quadruple  // 128-bit field data
};

#include <string>
#include <memory>

#include "error.hpp"

class FieldDescr {

  /// @class    FieldDescr
  /// @ingroup  Field
  /// @brief    Interface for the FieldDescr class

public: // public

  /// Initialize a FieldDescr object
  FieldDescr() throw();

  /// Destructor
  ~FieldDescr() throw();

  /// Copy constructor
  FieldDescr(const FieldDescr & field_descr) throw();

  /// Assignment operator
  FieldDescr & operator= (const FieldDescr & field_descr) throw();

  /// Return the number of fields
  int field_count() const throw();

  /// Return the integer handle for the named field
  int field_id(const std::string name) const throw();

  /// Return name of the ith field
  std::string field_name(size_t id_field) const throw()
  { return (id_field < field_name_.size()) ? field_name_[id_field] : ""; };

  /// Return the number of groups
  int group_count() const throw()
  { return group_name_.size(); };

  /// Return the integer handle for the named group
  int group_id(const std::string name) const throw();

  /// Return name of the ith group
  std::string group_name(int id_group) const throw();

  /// Return whether the given field is in the given group
  bool in_group(int id_field) const throw();

  /// alignment in bytes of fields in memory
  int alignment() const throw();

  /// padding in bytes between fields in memory
  int padding() const throw();

  /// centering of given field
  void centering(bool * cx, bool * cy, bool * cz) const throw();

  /// depth of ghost zones of given field
  void ghosts(int * gx, int * gy, int * gz) const throw();

  /// precision of given field
  precision_type precision() const throw();

  //----------------------------------------------------------------------

  // /// Return Field's name
  // std::string name () const throw();

  // /// Return centering of Field
  // const bool * centering () const throw();

  // /// Return centering of Field for the given axis
  // void set_centering (int axis, bool value) throw()
  // {
  //   if (0 <= axis && axis < 3) 
  //     centering_[axis] = value; 
  // };

  // /// Return Field minimum value
  // double min_value () const throw()
  // { return min_value_; };

  // /// Set Field minimum value
  // void set_min_value (double value) throw()
  // { min_value_ = value; };

  // /// Return Field maximum value
  // double max_value () const throw()
  // { return max_value_; };

  // /// Set Field maximum value
  // void set_max_value (double value) throw()
  // { max_value_ = value; };

  // /// Return action on violating Field minimum value
  // field_action min_action () const throw()
  // { return min_action_; };

  // /// Set action on violating Field minimum value
  // void set_min_action (field_action action) throw()
  // { min_action_ = action;  };

  // /// Return action on violating Field maximum value
  // field_action max_action () const throw()
  // { return max_action_; };

  // /// Set action on violating Field maximum value
  // void set_max_action (field_action action) throw()
  // { max_action_ = action;  };

  // /// Return precision of Field
  // precision_type precision () const throw()
  // { return precision_; };

  // /// Set precision of Field
  // void set_precision (precision_type precision) throw()
  // { 
  //   precision_ = (precision == precision_default) ? 
  //     default_precision_() : precision;
  // };

private: // functions

  precision_type default_precision_ () {

#if defined CONFIG_PRECISION_HALF
    return precision_half;
#elif defined CONFIG_PRECISION_SINGLE
    return precision_single;
#elif defined CONFIG_PRECISION_DOUBLE
    return precision_double;
#elif defined CONFIG_PRECISION_EXTENDED
    return precision_extended;
#elif defined CONFIG_PRECISION_QUADRUPLE
    return precision_quadruple;
#else
    return precision_unknown;
#endif

  };

private: // attributes

  /// Number of fields
  int field_count_;

  /// Number of field groups
  int group_count_;

  /// alignment of start of each field in bytes
  int alignment_;

  /// padding between fields in bytes
  int padding_;

  /// Courant number for fields
  double courant_;

  /// String defining each field
  std::vector<std::string> field_name_;

  /// String defining each group
  std::vector<std::string> group_name_;

  /// Precision of each field
  std::vector<precision_type> precision_;

  /// cell centering for each field
  std::vector<bool[3]> centering_;

  /// Ghost depth of each field
  std::vector<int[3]> ghosts_;

  /// minimum allowed value for each field
  std::vector<double> min_value_;

  /// maximum allowed value for each field
  std::vector<double> max_value_;

  /// what should be done if a field violates its minimum value
  std::vector<field_action> min_action_;

  /// what should be done if a field violates its maximum value
  std::vector<field_action> max_action_;
};

#endif /* FIELD_FIELD_DESCR_HPP */

