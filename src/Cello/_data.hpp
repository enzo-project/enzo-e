// See LICENSE_CELLO file for license and copyright information

/// @file     _data.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    Private include file for the \ref Data component

#ifndef _DATA_HPP
#define _DATA_HPP

enum particle_attribute_type {
  particle_attribute_unknown,
  particle_attribute_position,
  particle_attribute_velocity,
  particle_attribute_mass,
  particle_attribute_id,
  num_particle_attribute
};

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "pup_stl.h"

#include "_error.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

class FieldData;
class FieldFace;

#include "data_Grouping.hpp"
#include "data_ItField.hpp"
#include "data_ItFieldList.hpp"
#include "data_ItFieldRange.hpp"
#include "data_FieldDescr.hpp"
#include "data_FieldData.hpp"
#include "data_FieldFace.hpp"
#include "data_Field.hpp"

#include "data_ParticleDescr.hpp"
#include "data_ParticleBlock.hpp"

#include "data_Data.hpp"

#endif /* _DATA_HPP */

