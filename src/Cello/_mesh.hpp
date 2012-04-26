// See LICENSE_CELLO file for license and copyright information

/// @file     _mesh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 11 17:20:03 PST 2010
/// @brief    Private include file for the \ref Mesh component 

#ifndef _MESH_HPP
#define _MESH_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#ifdef CONFIG_USE_CHARM
#  include "pup_stl.h"
#endif

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

// Tree and components
#include "mesh_Node.hpp"
#include "mesh_NodeTrace.hpp"
#include "mesh_Tree.hpp"

// Hierarchy and components
#include "mesh_Block.hpp"
#include "mesh_Patch.hpp"
#include "mesh_Hierarchy.hpp"
#include "mesh_Factory.hpp"

// Iterators
#include "mesh_It.hpp"
#include "mesh_ItPatch.hpp"
#include "mesh_ItBlock.hpp"
#include "mesh_ItNode.hpp"

#endif /* _MESH_HPP */

