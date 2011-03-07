// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_HPP
#define MESH_HPP

/// @file     mesh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 11 17:20:03 PST 2010
/// @brief    Include file for the \ref Mesh component 

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "field.hpp"
#include "parallel.hpp"
#include "memory.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "mesh_Node2K.hpp"
#include "mesh_Node3K.hpp"
#include "mesh_TreeK.hpp"
#include "mesh_Tree2K.hpp"
#include "mesh_Tree3K.hpp"

#include "mesh_Block.hpp"

#include "mesh_Patch.hpp"
#include "mesh_Mesh.hpp"

#include "mesh_ItBlocks.hpp"


#endif /* MESH_HPP */

