// See LICENSE_CELLO file for license and copyright information

/// @file     _data.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    Private include file for the \ref Data component

#ifndef _DATA_HPP
#define _DATA_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

//----------------------------------------------------------------------
// Defines
//----------------------------------------------------------------------

#define PARTICLE_ALIGN 16

// integer limits on particle position within a Block:
//
//  -N    -N/2   0    N/2    N
//   *-----*-----*-----*-----*
//  LEFT   | IN BLOCK  | RIGHT
//        PMIN        PMAX
//   Full integer range is [-N,N)
//   Particles in block [-N/2,N/2)

#define PMIN_8 -64
#define PMAX_8  64

#define PMIN_16 -16384
#define PMAX_16  16384

#define PMIN_32 -1073741824
#define PMAX_32  1073741824

#define PMIN_64 -4611686018427387904
#define PMAX_64  4611686018427387904

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

class FieldData;
class FieldFace;

#include "data_Grouping.hpp"

#include "data_FieldDescr.hpp"
#include "data_FieldData.hpp"
#include "data_Field.hpp"
#include "data_FieldFace.hpp"

#include "data_ItIndex.hpp"
#include "data_ItIndexList.hpp"
#include "data_ItIndexRange.hpp"

#include "data_ParticleDescr.hpp"
#include "data_ParticleData.hpp"
#include "data_Particle.hpp"

#include "data_Data.hpp"

#include "data_DataMsg.hpp"

#endif /* _DATA_HPP */

