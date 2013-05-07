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

void CommBlock::p_adapt_enter()
{
  int initial_cycle     = simulation()->config()->initial_cycle;
  int initial_max_level = simulation()->config()->initial_max_level;

  count_adapt_ = (cycle() == initial_cycle) ? initial_max_level : 1;
  TRACE1("count_adapt = %d",count_adapt_);
  p_adapt(count_adapt_);
}

//----------------------------------------------------------------------

void CommBlock::p_adapt(int count)
{
  TRACE1("ADAPT p_adapt(%d)",count_adapt_);

  TRACE1("count_adapt = %d",count_adapt_);
  if (count_adapt_-- > 0) {

    TRACE1("count_adapt = %d",count_adapt_);
      int adapt = determine_adapt();

      TRACE1 ("ADAPT adapt = %s",adapt_name[adapt]);

      if      (adapt == adapt_refine)  p_refine();
      else if (adapt == adapt_coarsen) coarsen();

      CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt(), 
			    thisProxy[thisIndex]));

  } else {

    CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_exit(),
			  thisProxy[thisIndex]));

  }
}

//----------------------------------------------------------------------

int CommBlock::determine_adapt()
{
  TRACE("ADAPT CommBlock::determine_adapt()");

  // Only a leaf block can adapt

  if (! is_leaf()) return adapt_same;
  
  else {

    FieldDescr * field_descr = simulation()->field_descr();

    // Evaluate refinement criteria
    int i=0;
    int adapt = adapt_unknown;
    while (Refine * refine = simulation()->problem()->refine(i++)) {
      adapt = reduce_adapt_(adapt,refine->apply(this, field_descr));
    }

    // return the result
    return adapt;

  }

}

//----------------------------------------------------------------------

int CommBlock::reduce_adapt_(int a1, int a2) const throw()
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

void CommBlock::p_refine()
{
  TRACE("ADAPT CommBlock::p_refine()");

  const Factory * factory = simulation()->factory();
  int rank = simulation()->dimension();

  int nx,ny,nz;
  block()->field_block()->size(&nx,&ny,&nz);

  for (int ic=0; ic<num_child_; ic++) {

    int ic3[3];
    ic3[0] = (ic & 1) >> 0;
    ic3[1] = (ic & 2) >> 1;
    ic3[2] = (ic & 4) >> 2;

    TRACE5("ADAPT new child %d [%d %d %d]/ %d",
	   ic,ic3[0],ic3[1],ic3[2],num_child_);

    Index index = thisIndex;
    index.set_level(level_+1);
    index.set_child(level_+1,ic3[0],ic3[1],ic3[2]);
    index.clean();

    int num_field_blocks = 1;
    bool testing = false;

    // create new children
    factory->create_block 
      (thisProxy, index,
       nx,ny,nz,
       level_+1,
       num_field_blocks,
       count_adapt_,
       testing);

    // update own state
    child_exists_[ic] = true;

    // update neighbors' states, or invoke balance refine if not

    INCOMPLETE("CommBlock::p_refine(): "
	       "update neighbors' states, or invoke balance refine if not");
    for (int axis=0; axis<rank; axis++) {
      int face = ic3[axis];

      if (neighbor_exists_[IN(axis,face)]) {

	// if neighbor exists, tell it that it has a new nibling
	Index index = thisIndex;

      } else {

	int n3[3];
	size_forest (&n3[0],&n3[1],&n3[2]);
	Index index_uncle = index.index_uncle (axis,face,n3[axis]);
	thisProxy[index_uncle].p_refine();
	// if neighbor doesn't exist, then corresponding uncle it must refine
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
}

//----------------------------------------------------------------------

void CommBlock::q_adapt()
{
  thisProxy.doneInserting();
  TRACE1("count_adapt = %d",count_adapt_);
  TRACE3("ADAPT q_adapt %d %d %d",level_,count_adapt_,thisIndex.is_root());
  if (thisIndex.is_root()) {
    char buffer[40];
    TRACE1("count_adapt = %d",count_adapt_);
    // sprintf (buffer,"q_adapt(%d)",count_adapt_);
    // thisIndex.print(buffer);
    // thisProxy.p_print(buffer);
  thisProxy.p_adapt(count_adapt_);
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


