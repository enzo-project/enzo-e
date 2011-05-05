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

#include "cello.hpp"

#include "field.hpp"
#include "parallel.hpp"
#include "memory.hpp"
#include "simulation.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "mesh_tree.hpp"

class Factory;
#include "mesh_Block.hpp"

#include "mesh_Patch.hpp"
#include "mesh_Mesh.hpp"
#include "mesh_Factory.hpp"

#include "mesh_ItPatch.hpp"
#include "mesh_ItBlockLocal.hpp"


#endif /* MESH_HPP */

