// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/assorted/assorted.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for tests uncategorized source files within the \ref Enzo component
///
/// This subcomponent contains initializers for setting up test-problems

#ifndef ENZO_ASSORTED_ASSORTED_HPP
#define ENZO_ASSORTED_ASSORTED_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/mesh.hpp"    // Block
#include "Cello/problem.hpp" // Initial

#include "Enzo/enzo.hpp" // enzo_float

//----------------------------------------------------------------------
// Component headers
//----------------------------------------------------------------------

#include "assorted/EnzoMethodHeat.hpp"
#include "assorted/EnzoMethodM1Closure.hpp"
#include "assorted/EnzoMethodTurbulence.hpp"

#endif /* ENZO_ASSORTED_ASSORTED_HPP */
