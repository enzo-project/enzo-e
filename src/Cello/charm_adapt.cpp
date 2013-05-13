// See LICENSE_CELLO file for license and copyright information

/// @file     charm_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions

/// ADAPT
///
/// A mesh adaptation step involves evaluating refinement criteria
/// (Refine objects) on each leaf CommBlock in the hierarchy to
/// determine whether the CommBlock should refine, coarsen, or stay
/// the same.
///
/// REFINE
///
/// Each CommBlock tagged for refinement creates refined child
/// CommBlocks, and corresponding data are interpolated to the refined
/// CommBlocks.  Any CommBlock that is refined communicates this
/// information to its neighbors about its updated state.
///
/// COARSEN
///
/// Each CommBlock tagged for coarsening communicates with its parent,
/// which, if all children are tagged for coarsening, will coarsen by
/// deleting its children.  Corresponding data are restricted to
/// coarsened CommBlocks.  Any CommBlock that is coarsened 
/// tells its neighbors and parent about its updated state.
///
/// BALANCE
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

#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

const char * adapt_name[] = {
  "adapt_unknown",
  "adapt_same",
  "adapt_refine",
  "adapt_coarsen"
};

//======================================================================

void CommBlock::p_adapt_begin()
{
  Performance * performance = simulation()->performance();
  if (! performance->is_region_active(perf_adapt)) {
    performance->start_region(perf_adapt);
  }

  const Config * config = simulation()->config();

  int initial_cycle     = config->initial_cycle;
  int initial_max_level = config->initial_max_level;
  int mesh_max_level    = config->mesh_max_level;

  if (cycle() == initial_cycle) {
    count_adapt_ = initial_max_level;
  } else if (mesh_max_level > 0) {
    count_adapt_ = 1;
  } else {
    count_adapt_ = 0;
  }

  if (count_adapt_ > 0) {

    p_adapt(count_adapt_);

  } else {

    q_adapt_end();

  }
}

//----------------------------------------------------------------------

void CommBlock::p_adapt(int count)
{

  count_coarsen_ = 0;

  if (count_adapt_-- > 0) {

    int adapt = determine_adapt();

    if (adapt == adapt_refine)  p_refine();
    if (adapt == adapt_coarsen)   coarsen();

    CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt(), 
			  thisProxy[thisIndex]));

  } else {

    CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_end(),
			  thisProxy[thisIndex]));
  }
}

//----------------------------------------------------------------------

int CommBlock::determine_adapt()
{
  // Only leaves can adapt

  if (! is_leaf()) return adapt_same;

  else {

    FieldDescr * field_descr = simulation()->field_descr();

    int i=0;
    int adapt = adapt_unknown;

    Problem * problem = simulation()->problem();
    Refine * refine;

    while ((refine = problem->refine(i++))) {

      int adapt_step = refine->apply(this, field_descr);

      adapt = reduce_adapt_(adapt,adapt_step);

    }

    const Config * config = simulation()->config();
    int mesh_max_level    = config->mesh_max_level;

    if ((mesh_max_level == level_) && (adapt == adapt_refine)) {
      adapt = adapt_same;
    }

    return adapt;
  }
}

//----------------------------------------------------------------------

int CommBlock::reduce_adapt_(int a1, int a2) const throw()
{
  if (a1 == adapt_unknown) return a2;
  if (a2 == adapt_unknown) return a1;

  if ((a1 == adapt_coarsen) && 
      (a2 == adapt_coarsen)) {

    return adapt_coarsen;

  } else if ((a1 == adapt_refine)  || 
	     (a2 == adapt_refine)) {

    return adapt_refine;

  } else {

    return adapt_same;

  }
}

//----------------------------------------------------------------------

void CommBlock::p_refine()
{

  TRACE("CommBlock::p_refine()");
  int rank = simulation()->dimension();
  
  int nc = NC(rank);

  // block size
  int nx,ny,nz;
  block()->field_block()->size(&nx,&ny,&nz);

  // forest size
  int na3[3];
  size_forest (&na3[0],&na3[1],&na3[2]);

  int initial_cycle = simulation()->config()->initial_cycle;
  int balance       = simulation()->config()->mesh_balance;

  bool initial = initial_cycle == cycle();

  for (int ic=0; ic<nc; ic++) {

    int ic3[3];
    ic3[0] = (ic & 1) >> 0;
    ic3[1] = (ic & 2) >> 1;
    ic3[2] = (ic & 4) >> 2;

    Index index_child = index_.index_child(ic3);

    if ( ! is_child(index_child) ) {

      // create child

      int num_field_blocks = 1;
      bool testing = false;

      const Factory * factory = simulation()->factory();

      // create new children
      factory->create_block 
	(&thisProxy, index_child,
	 nx,ny,nz,
	 num_field_blocks,
	 count_adapt_,
	 initial,
	 testing);

      // update children
      set_child(index_child);

      int vc3[3];
      index_child.values(&vc3[0],&vc3[1],&vc3[2]);

      // for each adjacent neighbor
      for (int axis=0; axis<rank; axis++) {

	int face = ic3[axis];

	// update neighbor

	Index index_neighbor = index_.index_neighbor(axis,face,na3[axis]);

	if (is_neighbor(index_neighbor)) {

	  // new child is neighbor's nibling

	  thisProxy[index_neighbor].p_set_nibling(vc3);

	} else {

	  // no neighbor?  then uncle must refine

	  if (level_ > 0 && balance) {

	    Index index_uncle = index_.index_uncle (axis,face,na3[axis]);

	    thisProxy[index_uncle].p_refine();
	  }

	}

	// new child is nibling's neighbor

	Index index_nibling = index_.index_nibling
	  (axis, 1-face, ic3, na3[axis]);

	if (is_nibling(index_nibling)) {
	  thisProxy[index_nibling].p_set_neighbor(vc3);
	}

      }
    }
  }
  
}

//----------------------------------------------------------------------

void CommBlock::coarsen()
{
  if (level_ > 0) {
    Index index = thisIndex;
    int icx,icy,icz;
    index.child(level_,&icx,&icy,&icz);
    index.set_level(level_ - 1);
    index.clean();
    thisProxy[index].p_child_can_coarsen(IC(icx,icy,icz));
  }
}

//----------------------------------------------------------------------

void CommBlock::p_child_can_coarsen(int ic)
{
  INCOMPLETE1("CommBlock::p_child_can_coarsen(%d) is not implemented yet",
	      ic);
}

//----------------------------------------------------------------------

void CommBlock::q_adapt()
{
  thisProxy.doneInserting();

  if (thisIndex.is_root()) {
    thisProxy.p_adapt(count_adapt_);
  }
}

//----------------------------------------------------------------------

void CommBlock::q_adapt_end()
{
  TRACE("ADAPT CommBlock::q_adapt_end()");

  Performance * performance = simulation()->performance();
  if (performance->is_region_active(perf_adapt))
    simulation()->performance()->stop_region(perf_adapt);

  thisProxy.doneInserting();
  if (thisIndex.is_root()) {
    thisProxy.p_refresh_begin();
  };
  
}

//----------------------------------------------------------------------

void CommBlock::set_child(Index index)
{ children_.push_back(index); }

//----------------------------------------------------------------------
  
#endif /* CONFIG_USE_CHARM */


