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
///    CommBlock::p_refine()
///       parent.p_set_child_depth(my_level+1)
///       for neighbor
///          neighbor.p_set_neighbor_depth(index, my_level+1)
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
  TRACE1("ADAPT p_adapt(%d)",level);
  level_active_ = level;
  adapt (); 
}

//----------------------------------------------------------------------

void CommBlock::adapt()
{
  TRACE1("ADAPT adapt(%d)",level_active_);

  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();

  int max_level = simulation_charm->config()->mesh_max_level;

  if (max_level == 0) {

    refresh();

  } else {

    TRACE1("ADAPT CommBlock::adapt(%d)",level_active_);
    if (level_active_ == level_) {
      int adapt = determine_refine();
    }

    int max_level = simulation_charm->config()->mesh_max_level;
    if (level_active_ < max_level) {
      CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt(),thisProxy[thisIndex]));
    } else {
      refresh();
    }

  }

}

//----------------------------------------------------------------------

int CommBlock::determine_refine()
{
  TRACE("ADAPT CommBlock::determine_refine()");

  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();

  FieldBlock * field_block = block()->field_block();
  FieldDescr * field_descr = simulation_charm->field_descr();

  int i=0;
  int result = adapt_same;

  Refine * refine;
  while (refine = simulation_charm->problem()->refine(i++)) {
    result = refine->apply(field_block, field_descr);
    TRACE2("ADAPT refine %s %s",
	   refine->name().c_str(),
	   (result==adapt_refine) ? "refine" :
	   ((result==adapt_coarsen) ? "coarsen" : "same"));
  }
}

//----------------------------------------------------------------------

void CommBlock::p_refine()
{
  TRACE("ADAPT CommBlock::p_refine()");
}

//----------------------------------------------------------------------

void CommBlock::q_adapt()
{
  
  TRACE1("ADAPT q_adapt(%d)",level_active_);
  ++level_active_;
  adapt();
}

//----------------------------------------------------------------------

void CommBlock::p_coarsen()
{
  TRACE("ADAPT CommBlock::p_coarsen()");
}

//----------------------------------------------------------------------

void CommBlock::p_balance()
{
  TRACE("ADAPT CommBlock::p_balance()");
}

//----------------------------------------------------------------------


//======================================================================

#endif /* CONFIG_USE_CHARM */


