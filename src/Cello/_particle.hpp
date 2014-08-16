// See LICENSE_CELLO file for license and copyright information

/// @file     _particle.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    Private include file for the \ref Particle component

#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

/// @enum     particle_type
/// @brief    Mesh adaptation type: refine, coarsen, or stay the same

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

#include "pup_stl.h"

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "_error.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "field_Grouping.hpp"

#include "particle_ParticleDescr.hpp"
#include "particle_ParticleBlock.hpp"

#endif /* _PARTICLE_HPP */

