// $Id: field_FieldDescr.hpp 1261 2010-03-03 00:14:11Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELD_DESCR_HPP
#define FIELD_FIELD_DESCR_HPP

/// @file     field_FieldDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @todo     Replace ghosts/centering dynamic allocated arrays with vector to avoid big three
/// @todo     Split into "global" and "field-specific" attributes to reduce class size
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
  precision_quadruple  // 128-bit field data
};

#include <stdexcept>
#include <string>
#include <memory>
#include <set>

#include "error.hpp"

class FieldDescr {

  /// @class    FieldDescr
  /// @ingroup  Field
  /// @brief    Interface for the FieldDescr class

public: // functions

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

  /// Return name of the ith field
  std::string field_name(size_t id_field) const throw(std::out_of_range);

  /// Return the integer handle for the named field
  int field_id(const std::string name) const throw(std::out_of_range);

  /// Return the number of groups
  int group_count() const throw();

  /// Return name of the ith group
  std::string group_name(int id_group) const throw(std::out_of_range);

  /// Return the integer handle for the named group
  int group_id(const std::string name) const throw(std::out_of_range);

  /// Return whether the given field is in the given group
  bool field_in_group(int id_field, int id_group) const throw(std::out_of_range);

  /// alignment in bytes of fields in memory
  int alignment() const throw();

  /// padding in bytes between fields in memory
  int padding() const throw();

  /// courant number for fields
  double courant() const throw();

  /// centering of given field
  void centering(int id_field, bool * cx, bool * cy, bool * cz) const throw(std::out_of_range);

  /// depth of ghost zones of given field
  void ghosts(int id_field, int * gx, int * gy, int * gz) const throw(std::out_of_range);

  /// precision of given field
  precision_type precision(int id_field) const throw(std::out_of_range);

  /// Number of bytes per element required by the given field
  int bytes_per_element(int id_field) const throw();

  /// minimum value for the field
  double minimum_value(int id_field) const throw(std::out_of_range);

  /// minimum action for the field
  field_action minimum_action(int id_field) const throw(std::out_of_range);

  /// maximum value for the field
  double maximum_value(int id_field) const throw(std::out_of_range);

  /// maximum action for the field
  field_action maximum_action(int id_field) const throw(std::out_of_range);

  //----------------------------------------------------------------------

  /// Insert a new field
  void insert_field(std::string name_field) throw();

  /// Insert a new group
  void insert_group(std::string name_group) throw();

  /// Set membership of a field in a group
  void set_field_in_group(int id_field, int id_group) throw(std::out_of_range);

  /// Set alignment
  void set_alignment(int alignment) throw();

  /// Set padding
  void set_padding(int padding) throw();

  /// Set courant
  void set_courant(double courant) throw();

  /// Set centering for a field
  void set_centering(int id_field, bool cx, bool cy, bool cz) throw(std::out_of_range);

  /// Set ghosts for a field
  void set_ghosts(int id_field, int gx, int gy, int gz) throw(std::out_of_range);

  /// Set precision for a field
  void set_precision(int id_field, precision_type precision) throw(std::out_of_range);

  /// Set minimum bound and action
  void set_minimum (int id_field, double min_value, field_action min_action) 
    throw(std::out_of_range);

  /// Set maximum bound and action
  void set_maximum (int id_field, double max_value, field_action max_action) 
    throw(std::out_of_range);

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

  /// alignment of start of each field in bytes
  int alignment_;

  /// padding between fields in bytes
  int padding_;

  /// Courant number for fields
  double courant_;

  /// String defining each field
  std::vector<std::string> field_name_;

  /// Integer id for each field.  Inverse mapping of field_name_
  std::map<std::string,int> field_id_;

  /// String defining each group
  std::vector<std::string> group_name_;

  /// Integer id for each group.  Inverse mapping of group_name_
  std::map<std::string,int> group_id_;

  typedef std::set<int> set_int_type;
  /// Set of groups containing each field.  field_in_group_[field][group]
  std::vector<set_int_type> field_in_group_;

  /// Precision of each field
  std::vector<precision_type> precision_;

  /// cell centering for each field
  std::vector<bool *> centering_;

  /// Ghost depth of each field
  std::vector<int *> ghosts_;

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

