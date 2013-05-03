// See LICENSE_CELLO file for license and copyright information

/// @file     charm_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions
///
/// ADAPT

/// A mesh adaptation step involves evaluating refinement criteria
/// (Refine objects) on each leaf CommBlock in the hierarchy to
/// determine whether the CommBlock should refine, coarsen, or stay
/// the same.
///
/// Each CommBlock tagged for refinement creates refined child
/// CommBlocks.  Each CommBlock tagged for coarsening communicates
/// with its parent, which, if all children are tagged for coarsening,
/// will coarsen by deleting its children.  Corresponding data are
/// interpolated to refined CommBlocks and coarsened to coarsened
/// CommBlocks.
///
/// Any CommBlock that is coarsened or refined tells its neighbors and
/// parent about its updated state (number of descendents).  The
/// "depth" of each child (distance from deepest descendent) is stored
/// in each CommBlock, which determines the depth of the CommBlock.
///
/// After "quiescence" (wait until no communication), a balancing
/// phase is performed.  The mesh is traversed by levels, finest
/// first, and any CommBlock that is adjacent to any CommBlock that
/// has a grandchild is tagged for refinement.  Balancing a CommBlock
/// can trigger further CommBlocks to require balancing, but only
/// coarser ones, which will be handled in the next level.  A
/// quiescence step is used between each level.
///
/// Typically only one mesh adaptation step is performed at a time,
/// except when applying initial conditions.  In that case, several
/// steps may be applied, up to a specified maximum
/// (Mesh:initial_max_level).

   
/// CommBlock::p_adapt(level) 
///    if (adapt)
///       adapt(level)
///    else
///       refresh()
///
/// CommBlock::adapt(level)
///    if (level == my_level && is_leaf())
///       adapt = determine_refinement()
///       if (adapt == adapt_refine)  p_refine()
///       if (adapt == adapt_coarsen) p_coarsen()
///    if (level < max_level)
///       StartQD( p_adapt(level + 1) )
///    else
///       >>>>> refresh() >>>>>
///
/// CommBlock::p_refine()
///    parent.p_set_child_depth(my_level+1)
///    for neighbor
///       neighbor.p_set_neighbor_depth(index, my_level+1)
///       neighbor.p_balance(index)
///    for child
///       create child
///
/// CommBlock::p_coarsen()
///    ???
///
/// CommBlock::p_balance(index)
///    if (is_leaf())
///       adapt = determine_balance()
///       if (adapt == adapt_refine) p_refine()
///
/// CommBlock::determine_balance()
///    adapt_type = adapt_same     
///    for neighbor
///       if (neighbor_depth > my_level) 
///          adapt_type = adapt_refine
///    return adapt_type


#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// NOTE: +2 is so indices can be -1, e.g. for child indicies of neighboring blocks
#define IC(icx,icy,icz)  ((icx+2)%2) + 2*( ((icy+2)%2) + 2*((icz+2)%2) )
#define NC(rank) (rank == 1 ? 2 : (rank == 2 ? 4 : 8))
#define IN(axis,face)  ((face) + 2*(axis))
#define NN(rank) (rank*2)

const char * adapt_name[] = {
  "adapt_unknown",
  "adapt_same",
  "adapt_refine",
  "adapt_coarsen"
};

//======================================================================

void CommBlock::p_phase_adapt() {  p_adapt(0); }

//======================================================================

void CommBlock::p_adapt(int level)
{
  level_active_ = level;
  TRACE1("ADAPT p_adapt(%d)",level);
  adapt (); 
}

//----------------------------------------------------------------------

void CommBlock::adapt()
{
  TRACE1("ADAPT adapt(%d)",level_active_);

  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();

  int mesh_max_level    = simulation_charm->config()->mesh_max_level;
  int initial_max_level = simulation_charm->config()->initial_max_level;
  int initial_cycle     = simulation_charm->config()->initial_cycle;

  TRACE2 ("ADAPT cycle %d  initial_cycle %d",cycle(),initial_cycle);
  int max_level = (cycle() == initial_cycle) ? 
    initial_max_level : mesh_max_level;

  TRACE3("ADAPT max_level %d  initial %d  mesh %d",max_level,
	 initial_max_level, mesh_max_level);

  if (max_level == 0) {

    CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_exit(),thisProxy[thisIndex]));

  } else {

    TRACE2("ADAPT determine_refine %d %d",level_active_,level_);
    if (level_active_ == level_) {
      int adapt = determine_refine();
      TRACE1 ("ADAPT adapt = %s",adapt_name[adapt]);
      if      (adapt == adapt_refine)  refine();
      else if (adapt == adapt_coarsen) coarsen();
    }

    int max_level = simulation_charm->config()->mesh_max_level;

    if (level_active_ < max_level) {

      CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt(),thisProxy[thisIndex]));

    } else {

      level_active_ = -1;

      CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_exit(),thisProxy[thisIndex]));

    }
  }
}

//----------------------------------------------------------------------

int CommBlock::determine_refine()
{
  TRACE("ADAPT CommBlock::determine_refine()");

  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation_charm->field_descr();

  int i=0;
  int adapt = adapt_unknown;
  while (Refine * refine = simulation_charm->problem()->refine(i++)) {
    adapt = update_adapt_(adapt,refine->apply(this, field_descr));
  }
  return adapt;

}

//----------------------------------------------------------------------

int CommBlock::update_adapt_(int a1, int a2) const throw()
{
  TRACE2("ADAPT %d %d",a1,a2);
  if (a1 == adapt_unknown) return a2;
  if (a2 == adapt_unknown) return a1;
  if      ((a1 == adapt_coarsen) && (a2 == adapt_coarsen)) 
    return adapt_coarsen;
  else if ((a1 == adapt_refine)  || (a2 == adapt_refine))
    return adapt_refine;
  else
    return adapt_same;
}

//----------------------------------------------------------------------

void CommBlock::p_update_depth (int ic, int depth)
{
  // update self
  depth_[ic] = std::max(depth_[ic],depth);
  // update parent
  if (level_ > 0) {
    Index index = thisIndex;
    int icx,icy,icz;
    index.child(level_,&icx,&icy,&icz);
    index.set_level(level_-1);
    index.clean();
    thisProxy[index].p_update_depth (IC(icx,icy,icz), depth + 1);
  }
}

//----------------------------------------------------------------------

void CommBlock::refine()
{
  TRACE("ADAPT CommBlock::refine()");

  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();
  const Factory * factory = simulation_charm->factory();
  int rank = simulation_charm->dimension();

  int ibx = index_[0];
  int iby = index_[1];
  int ibz = index_[2];
  TRACE3("ADAPT ib = %d %d %d",ibx,iby,ibz);
  int nbx = size_[0];
  int nby = size_[1];
  int nbz = size_[2];

  int nx,ny,nz;
  block()->field_block()->size(&nx,&ny,&nz);

  double xm[2],ym[2],zm[2];
  block()->lower(&xm[0],&ym[0],&zm[0]);
  double xp[2],yp[2],zp[2];
  block()->upper(&xp[1],&yp[1],&zp[1]);

  xm[1] = xp[0] = 0.5 * (xm[0] + xp[1]);
  ym[1] = yp[0] = 0.5 * (ym[0] + yp[1]);
  zm[1] = zp[0] = 0.5 * (zm[0] + zp[1]);

  if (! ((0.0 <= xm[0] && xm[0] <= 1.0) &&
	 (0.0 <= ym[0] && ym[0] <= 1.0) &&
	 (0.0 <= zm[0] && zm[0] <= 1.0) &&
	 (0.0 <= xm[1] && xm[1] <= 1.0) &&
	 (0.0 <= ym[1] && ym[1] <= 1.0) &&
	 (0.0 <= zm[1] && zm[1] <= 1.0)) ) {
    TRACE6("CommBlock ERROR:  %f %f %f  %f %f %f", xm[0],ym[0],zm[0],xp[0],yp[0],zp[0]);
    TRACE6("CommBlock ERROR:  %f %f %f  %f %f %f", xm[1],ym[1],zm[1],xp[1],yp[1],zp[1]);
  }
  int num_field_blocks = 1;
  bool testing = false;

  int nc = NC(rank);

  for (int ic=0; ic<nc; ic++) {

    int icx = (ic & 1) >> 0;
    int icy = (ic & 2) >> 1;
    int icz = (ic & 4) >> 2;

    TRACE5("ADAPT new child %d [%d %d %d]/ %d",ic,icx,icy,icz,nc);

    Index index = thisIndex;
    index.set_level(level_+1);
    index.set_tree (level_+1,icx,icy,icz);
    index.clean();

    factory->create_block 
      (thisProxy, index,
       nbx,nby,nbz,
       nx,ny,nz,
       level_+1,
       xm[icx],ym[icy],zm[icz],
       xp[icx],yp[icy],zp[icz],
       num_field_blocks,
       testing);

    p_update_depth (ic,1);

  }
  
}

//----------------------------------------------------------------------

void CommBlock::coarsen()
{
  TRACE("ADAPT CommBlock::coarsen()");
}

//----------------------------------------------------------------------

void CommBlock::q_adapt()
{
  thisProxy.doneInserting();
  TRACE3("ADAPT q_adapt %d %d %d",level_,level_active_,thisIndex.is_root());
  if (thisIndex.is_root()) {
    char buffer[40];
    sprintf (buffer,"q_adapt(%d)",level_active_);
    thisIndex.print(buffer);
    thisProxy.p_print(buffer);
    thisProxy.p_adapt(level_active_ + 1);
  }
}

//----------------------------------------------------------------------

void CommBlock::p_balance()
{
  TRACE("ADAPT CommBlock::p_balance()");
}

//----------------------------------------------------------------------

void CommBlock::q_adapt_exit()
{
  TRACE("ADAPT CommBlock::q_adapt_exit()");
  thisProxy.doneInserting();
  if (thisIndex.is_root()) {
    thisProxy.p_refresh();
  };
  
}

#endif /* CONFIG_USE_CHARM */


