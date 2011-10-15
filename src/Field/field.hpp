// See LICENSE_CELLO file for license and copyright information

/// @file     field.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    Include file for the \ref Field component

#ifndef FIELD_HPP
#define FIELD_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "pngwriter.h"

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "error.hpp"
#include "disk.hpp"
#include "mesh.hpp"
// #include "method.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

class FieldBlock;
class FieldFace;

#include "field_ItField.hpp"
#include "field_ItFieldList.hpp"
#include "field_ItFieldRange.hpp"
#include "field_FieldDescr.hpp"
#include "field_FieldBlock.hpp"
#include "field_FieldFace.hpp"
#include "field_FacesGroup.hpp"

#endif /* FIELD_HPP */

