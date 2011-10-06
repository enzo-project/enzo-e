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
  :  count_output_(0),
     count_prepare_(0),
     dt_hierarchy_(std::numeric_limits<double>::max()),
     stop_hierarchy_(true)
{
  TRACE("BlockReduce()");
}

#endif

//----------------------------------------------------------------------

// See Simulation/simulation_charm_output.cpp for BlockReduce::p_prepare()

#ifdef CONFIG_USE_CHARM

void BlockReduce::p_prepare(int    count, 
			    int    cycle, 
			    double time,
			    double dt_block, 
			    int    stop_block)
{

  TRACE1("BlockReduce::p_prepare(%d)",count);
  // Parallel reduction of dt and stopping criteria 

  dt_hierarchy_   = MIN(dt_hierarchy_, dt_block);
  stop_hierarchy_ = stop_hierarchy_ && stop_block;

  // Assumes cycle and time are the same for all "incoming" blocks;
  // only use the last one

  if (++count_prepare_ >= count) {

    //--------------------------------------------------
    // Simulation::p_output()
    //--------------------------------------------------

    TRACE1("BlockReduce::p_prepare(%d) calling Simulation::p_output()",count);
    proxy_simulation.p_output(cycle, time, dt_hierarchy_, stop_hierarchy_);

    // Reset pool

    count_prepare_ = 0;
    dt_hierarchy_ = std::numeric_limits<double>::max();
    stop_hierarchy_ = true;

  }
}

#endif

//----------------------------------------------------------------------

// SEE simulation_charm_output.cpp for
// BlockReduce::p_output_reduce()

//----------------------------------------------------------------------
