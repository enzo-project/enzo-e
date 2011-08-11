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
     dt_mesh_(std::numeric_limits<double>::max()),
     stop_mesh_(true)
{
}

//----------------------------------------------------------------------

void BlockReduce::p_prepare(int    count, 
				 int    cycle, 
				 double time,
				 double dt_block, 
				 int    stop_block)
{
  TRACE("BlockReduce::p_prepare");
  // Assumes cycle and time are the same for all "incoming" blocks;
  // only use the last one

  //--------------------------------------------------
  // Timestep
  //--------------------------------------------------

  dt_mesh_   = MIN(dt_mesh_, dt_block);

  //--------------------------------------------------
  // Stopping
  //--------------------------------------------------

  stop_mesh_ = stop_mesh_ && stop_block;

  if (++count_prepare_ >= count) {

    //--------------------------------------------------
    // Simulation::p_output()
    //--------------------------------------------------
    proxy_simulation.p_output(cycle, time, dt_mesh_, stop_mesh_);

    // Reset pool
    count_prepare_ = 0;
    dt_mesh_ = std::numeric_limits<double>::max();
    stop_mesh_ = true;

  }
}

#endif
//----------------------------------------------------------------------

//  --- Accumulate block output contributions and write output to disk ---

#ifdef CONFIG_USE_CHARM

void BlockReduce::p_output_reduce(int count)
{
  TRACE("BlockReduce::p_output_reduce");
  if (++count_output_ >= count) {
    INCOMPLETE("BlockReduce::p_output_reduce()");
    proxy_simulation.p_output_reduce();
    count_output_ = 0;
  }
}


#endif
