// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineRegion.cpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     2019-04-13
/// @brief    [\ref Mesh] Implementation of the RefineRegion class
///

#include "mesh.hpp"

//---------------------------------------------------------------------

RefineRegion::RefineRegion() throw()
{
  TRACE("RefineRegion::RefineRegion");
}

//----------------------------------------------------------------------

int RefineRegion::apply( Block * block ) throw()
{
  // This does nothing. RefineRegion is enforced in control_adapt
  // since this refinement will override any other refinement criteria
  // used. This allows a check to see if a given block is in a static
  // refinement region to see if any other refinement criteria need to be
  // checked.

  int adapt_result = 0;
  return adapt_result;
}

//=========================================================================
