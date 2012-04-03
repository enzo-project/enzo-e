// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_BlockReduce.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-08-10
/// @brief    Implementation of the BlockReduce class

#include "cello.hpp"

#include "mesh.hpp"
#include "simulation.hpp"

#include "mesh_charm.hpp"
#include "simulation_charm.hpp"
//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

BlockReduce::BlockReduce()
  : count_output_(0)
{
}

#endif

//----------------------------------------------------------------------

// See simulation_charm_output.cpp for BlockReduce::entry_output()

// See simulation_charm_output.cpp for BlockReduce::entry_output_reduce()
