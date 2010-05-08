// $Id: field_FieldDescr.hpp 1261 2010-03-03 00:14:11Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELD_DESCR_HPP
#define FIELD_FIELD_DESCR_HPP

/// @file     field_FieldDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    Declaration for the FieldDescr class

enum action_field {
  action_field_unknown,  // Uninitialized action
  action_field_none,     // Do nothing if range exceeded
  action_field_set,      // Set field values to min / max if range exceeded
  action_field_warning,  // Issue warning if range exceeded
  action_field_error,    // Issue error if range exceeded
  action_field_timestep, // Retry with reduced timestep if range exceeded
  action_field_method    // Retry with alternate method if range exceeded
};


enum precision_type {
  precision_unknown,  //  unknown precision
  precision_32bit,    //  32-bit field data
  precision_64bit,    //  64-bit field data
};

#include <string>
#include <memory>

class FieldDescr {

  /// @class    FieldDescr
  /// @ingroup  Field
  /// @brief    Interface for the FieldDescr class

public: // public

  /// Initialize a FieldDescr object
  FieldDescr();

  /// Set Field name
  void set_name (std::string name)
  { name_ = name; };

  /// Get Field name
  std::string name ()
  { return name_; };

private: // attributes

  /// String defining the field's name
  std::string name_;

  /// Integer handle identifying the Field (somewhat redundant with
  /// block_* attributes)
  int id_;

  /// Dimension of the Field
  int dim_;

  /// Identify which Block contains the Field data
  int block_number_;

  /// Identify where in the Block is the Field data
  int block_offset_ ;

  /// Cell centering, defined as (0,0,0) <= (px,py,pz) <= (1,1,1)
  /// Cell centered = (.5,.5,.5)
  std::auto_ptr<float> centering_;

  /// Minimum allowed value for the Field
  double min_;

  /// Maximum allowed value for the Field
  double max_;

  /// Action to perform if Field values go below min_
  action_field min_action_ ;

  /// Action to perform if Field values go above max_
  action_field max_action_ ;

  /// Precision of the Field data
  precision_type precision_;

};

#endif /* FIELD_FIELD_DESCR_HPP */

