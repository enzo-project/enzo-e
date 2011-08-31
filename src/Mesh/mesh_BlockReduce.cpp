// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_BlockReduce.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file mesh_BlockReduce.cpp
///
/// Detailed description of file mesh_BlockReduce.cpp

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
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void BlockReduce::p_prepare(int    count, 
			    int    cycle, 
			    double time,
			    double dt_block, 
			    int    stop_block)
{
  // Assumes cycle and time are the same for all "incoming" blocks;
  // only use the last one

  //--------------------------------------------------
  // Timestep
  //--------------------------------------------------

  dt_hierarchy_   = MIN(dt_hierarchy_, dt_block);

  //--------------------------------------------------
  // Stopping
  //--------------------------------------------------

  stop_hierarchy_ = stop_hierarchy_ && stop_block;

  if (++count_prepare_ >= count) {

    //--------------------------------------------------
    // Simulation::p_output()
    //--------------------------------------------------
    proxy_simulation.p_output(cycle, time, dt_hierarchy_, stop_hierarchy_);

    // Reset pool
    count_prepare_ = 0;
    dt_hierarchy_ = std::numeric_limits<double>::max();
    stop_hierarchy_ = true;

  }
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void BlockReduce::p_output_reduce(int count)
{
  if (++count_output_ >= count) {
    INCOMPLETE("BlockReduce::p_output_reduce()");
    proxy_simulation.p_output_reduce();
    count_output_ = 0;
  }
}

#endif

//----------------------------------------------------------------------
