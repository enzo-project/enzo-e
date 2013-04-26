// See LICENSE_CELLO file for license and copyright information

/// @file     charm_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions
///
///    ADAPT
///
///    CommBlock::p_adapt(level) 
///       if (adapt)
///          adapt(level)
///       else
///          refresh()
///
///    CommBlock::adapt(level)
///       if (level == my_level && is_leaf())
///          adapt = determine_refinement()
///          if (adapt == adapt_refine)  p_refine()
///          if (adapt == adapt_coarsen) p_coarsen()
///       if (level < max_level)
///          StartQD( p_adapt(level + 1) )
///       else
///          >>>>> refresh() >>>>>
///
///    CommBlock::p_refine(level)
///       parent.p_set_child_depth(level+1)
///       for neighbor
///          neighbor.p_set_neighbor_depth(index, level+1)
///          neighbor.p_balance(index)
///       for child
///          create child
///
///    CommBlock::p_coarsen()
///       ???
///
///    CommBlock::p_balance(index)
///       if (is_leaf())
///          adapt = determine_balance()
///          if (adapt == adapt_refine) p_refine()
///
///    CommBlock::determine_balance()
///       adapt_type = adapt_same     
///       for neighbor
///          if (neighbor_depth > my_level) 
///             adapt_type = adapt_refine
///       return adapt_type


#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void CommBlock::p_adapt(int level)
{
  TRACE("p_adapt()");
  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();

  if (simulation_charm->config()->mesh_max_level == 0) {

    refresh();

  } else {

    INCOMPLETE("CommBlock::p_adapt(): Adaptation not implemented, calling refresh()");

    refresh();

  }
  
}
//======================================================================

#endif /* CONFIG_USE_CHARM */


