// See LICENSE_CELLO file for license and copyright information

/// @file     _mesh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 11 17:20:03 PST 2010
/// @brief    Private include file for the \ref Mesh component 

#ifndef _MESH_HPP
#define _MESH_HPP

/// @enum     adapt_type
/// @brief    Mesh adaptation type: refine, coarsen, or stay the same

enum adapt_type {
  adapt_unknown,
  adapt_coarsen,
  adapt_same,
  adapt_refine
};

enum refresh_type {
  refresh_unknown,
  refresh_coarse,
  refresh_same,
  refresh_fine
};

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include "pup_stl.h"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

// Hierarchy and components

class Tree;

#include "mesh_Index.hpp"

#include "mesh_Block.hpp"
#include "mesh_Hierarchy.hpp"
#include "mesh_Factory.hpp"

// Tree and components
#include "mesh_Node.hpp"
#include "mesh_NodeTrace.hpp"
#include "mesh_Tree.hpp"

// Iterators
#include "mesh_It.hpp"
#include "mesh_ItBlock.hpp"
#include "mesh_ItNode.hpp"

// Refinement
#include "mesh_Refine.hpp"
#include "mesh_RefineSlope.hpp"
#include "mesh_RefineMass.hpp"
#include "mesh_RefineMask.hpp"
#include "mesh_ItFace.hpp"
#include "mesh_ItChild.hpp"

#endif /* _MESH_HPP */

