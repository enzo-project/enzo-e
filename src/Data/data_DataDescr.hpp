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
  DataDescr() throw();

  //----------------------------------------------------------------------
  // Field functions
  //----------------------------------------------------------------------

  /// Return the Field descriptor

  FieldDescr * field_descr ();

  //----------------------------------------------------------------------
  // Particle functions
  //----------------------------------------------------------------------

private: // attributes

  int           dim_;

  FieldDescr    field_descr_;
  ParticleDescr particle_descr_;

};

#endif /* DATA_DATA_HPP */

