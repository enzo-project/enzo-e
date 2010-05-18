// $Id: data_data.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    Declaration of the DataDescr class

#ifndef DATA_DATA_HPP
#define DATA_DATA_HPP

#include <vector>

#include "particle.hpp"
#include "field.hpp"

class DataDescr {

  /// @class    DataDescr
  /// @ingroup  Data
  /// @brief    Container class for Particle and Field descriptors

public: // interface

  /// Initialize the DataDescr object
  DataDescr()
    : dim_(3),
      field_descr_(),
      particle_descr_()
  { };


  /// Set dimensions
  void set_dimension (int dim)
  {
    if (1 <= dim && dim <= 3) {
      dim_ = dim;
    } else {
      ERROR_MESSAGE("DataDescr::set_dimension","dim out of range");
    }
  }

  int dimension () const throw ()
  { return dim_; }

  //----------------------------------------------------------------------
  // Field functions
  //----------------------------------------------------------------------

  /// Return the ith Field descriptor
  FieldDescr * field_descr (int index)
  { 
    return (0 <= index && index < field_count()) ? 
      field_descr_[index] : 0;
  };

  /// Add a new field to the list of known fields
  int field_insert (std::string name) throw()
  { 
    int index = field_count();
    FieldDescr * field = new FieldDescr(name,dim_);
    field_descr_.push_back(field); 
    return index;
  };

  /// Return the number of Fields
  int field_count() const throw()
  { return field_descr_.size(); }


  /// Return named Field's index
  int field_index (std::string name) const throw()
  { 
    for (int index=0; index < field_count(); index++) {
      if (field_descr_[index]->name() == name) return index;
    }
    return -1; // Uh oh
  }

  /// Return indexed Field's name
  std::string field_name (int index) const throw()
  { 
    return (0 <= index && index < field_count()) ? 
      field_descr_[index]->name() : "";
  }

  //----------------------------------------------------------------------
  // Particle functions
  //----------------------------------------------------------------------

private: // attributes

  int           dim_;

  std::vector<FieldDescr *>    field_descr_;
  std::vector<ParticleDescr *> particle_descr_;

};

#endif /* DATA_DATA_HPP */

