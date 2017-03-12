// See LICENSE_CELLO file for license and copyright information

/// @file     charm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-20 16:37:14
/// @brief    Include file for the \ref Charm component 

#include "charm++.h"

#include <string>

#include "_error.hpp"
#include "mesh_Index.hpp"
#include "charm_reductions.hpp"
#include "charm_Call.hpp"
#include "charm_Sync.hpp"
#include "charm_MappingArray.hpp"
#include "charm_MappingTree.hpp"

#include "charm_MsgRefresh.hpp"
#include "charm_MsgCoarsen.hpp"
#include "charm_MsgRefine.hpp"
