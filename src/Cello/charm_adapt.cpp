// See LICENSE_CELLO file for license and copyright information

/// @file     charm_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions

/// Maintain type (coarse,same,fine) of neighbor at all required faces
/// Update when refining
/// Update when coarsening
/// Force balanced-refined shouldn't coarsen
///     determine adapt
///     neighbor update
///     balance check
///
/// adapt: refine / coarsen at all levels
///
/// 1. determine adapt
/// 2. refine / balance where necessary
/// 3. QD
/// 4. coarsen where possible
/// 5. QD [optional?]
///
/// p_adapt_begin()
///    adapt_ = determine_adapt 
///    if (refine_) {
///      tell known neighbors
///      create children
///      update children neighbors
///      balance if necessary
///    }
///  each face stores relative level of finest grid
///  0==same, -1==coarse, +1== fine


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
  TRACE ("BEGIN PHASE ADAPT");
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

    p_adapt_start();

  } else {

    q_adapt_end();

  }
}

//----------------------------------------------------------------------

void CommBlock::p_adapt_start()
{
  // Initialize child coarsening counter
  count_coarsen_ = 0;

  // Inhibit coarsening on startup
  int initial_cycle = simulation()->config()->initial_cycle;

  if (count_adapt_-- > 0) {

    adapt_ = determine_adapt();

    if (adapt_ == adapt_refine)  refine();

    CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_next(), 
			  thisProxy[thisIndex]));

  } else {

    CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_end(),
			  thisProxy[thisIndex]));
  }
}

//----------------------------------------------------------------------

void CommBlock::q_adapt_next()
{
  if (adapt_ == adapt_coarsen) coarsen();

  CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_stop(), 
			thisProxy[thisIndex]));
}

//----------------------------------------------------------------------

int CommBlock::determine_adapt()
{
  // Only leaves can adapt

  if (! is_leaf()) return adapt_same;

  else {

    FieldDescr * field_descr = simulation()->field_descr();

    int index_refine = 0;
    int adapt = adapt_unknown;

    Problem * problem = simulation()->problem();
    Refine * refine;

    while ((refine = problem->refine(index_refine++))) {

      int adapt_step = refine->apply(this, field_descr);

      adapt = reduce_adapt_(adapt,adapt_step);

    }

    const Config * config = simulation()->config();
    int mesh_max_level    = config->mesh_max_level;

    if ((adapt == adapt_refine)  && ! can_refine())  adapt = adapt_same;
    if ((adapt == adapt_coarsen) && ! can_coarsen()) adapt = adapt_same;

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

void CommBlock::refine()
{

  adapt_ = adapt_unknown;

#ifdef CELLO_TRACE
  index_.print ("refine",-1,2);
#endif /* CELLO_TRACE */
  
  TRACE("CommBlock::refine()");
  int rank = simulation()->dimension();
  
  int nc = NC(rank);

  // block size
  int nx,ny,nz;
  block()->field_block()->size(&nx,&ny,&nz);

  // forest size
  int na3[3];
  size_forest (&na3[0],&na3[1],&na3[2]);

  int initial_cycle = simulation()->config()->initial_cycle;
  bool do_balance   = simulation()->config()->mesh_balance;

  bool initial = (initial_cycle == cycle());

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

      // <duplicated code: refactor me!>
      int narray = 0;  
      char * array = 0;
      int op_array = op_array_prolong;

      Simulation * simulation = proxy_simulation.ckLocalBranch();
      FieldDescr * field_descr = simulation->field_descr();
      FieldBlock * field_block = block_->field_block();
      Prolong * prolong = simulation->problem()->prolong();

      FieldFace field_face (field_block,field_descr);

      field_face.set_prolong(prolong,ic3[0],ic3[1],ic3[2]);
      field_face.set_face(0,0,0); // 0-face is full block

      field_face.load(&narray, &array);

      // </duplicated code>
    
      const Factory * factory = simulation->factory();

      factory->create_block 
	(&thisProxy, index_child,
	 nx,ny,nz,
	 num_field_blocks,
	 count_adapt_,
	 initial,
	 cycle_,time_,dt_,
	 narray, array, op_array,
	 testing);

      set_child(index_child);

      // for each (axis,face)

      for (int axis=0; axis<rank; axis++) {

	int face = ic3[axis];

	int in3[3] = {0};
	in3[axis] = 2*face-1; // (0,1) --> (+1,-1)

	// update non-sibling neighbors

	Index index_neighbor = index_.index_neighbor(axis,face,na3[axis]);

	if (is_neighbor(in3)) {

	  // new child is neighbor's nibling

	  thisProxy[index_neighbor].p_set_nibling (in3);

	} else {

	  // no neighbor?  then uncle must refine

	  if (level_ > 0 && do_balance) {

	    Index index_uncle = index_neighbor.index_parent();

	    thisProxy[index_uncle].p_balance();
	  }

	}

	// new child is nibling's neighbor (and vice versa)

	int jc3[3] = {ic3[0],ic3[1],ic3[2]};
	jc3[axis] = 1-jc3[axis];

	if (is_nibling(ic3)) {
	  
	  Index index_nibling = index_neighbor.index_child(jc3);

	  int in3[3] = {0};
 	  in3[axis] = 2*face - 1; // (0,1) --> (+1,-1)

	  thisProxy[index_child].p_set_neighbor (in3);

	  in3[axis] = -in3[axis];

	  thisProxy[index_nibling].p_set_neighbor (in3);
	}

      }
    }
  }
  
}

//----------------------------------------------------------------------

void CommBlock::p_balance()
{
  refine();
}

//----------------------------------------------------------------------

bool CommBlock::can_refine() const
{
  int max_level = simulation()->config()->mesh_max_level;
  return (level_ < max_level);
}

//----------------------------------------------------------------------

bool CommBlock::can_coarsen() const
{ 

  if (level_ <= 0) return false;

  int refresh_rank = simulation()->config()->field_refresh_rank;
  int rank = simulation()->dimension();
  int if3m[3],if3p[3],if3[3];

  loop_limits_refresh_(if3m+0,if3m+1,if3m+2,
		       if3p+0,if3p+1,if3p+2);

  for (if3[0]=if3m[0]; if3[0]<=if3p[0]; if3[0]++) {
    for (if3[1]=if3m[1]; if3[1]<=if3p[1]; if3[1]++) {
      for (if3[2]=if3m[2]; if3[2]<=if3p[2]; if3[2]++) {
	int face_rank = rank - (abs(if3[0]) + abs(if3[1]) + abs(if3[2]));
	if (refresh_rank <= face_rank && face_rank < rank) {
	  int i=IN3(if3);
	  if (face_level_[i] > level_) return false;
	}
      }
    }
  }

  return true;
}

//----------------------------------------------------------------------

void CommBlock::coarsen()
{

  adapt_ = adapt_unknown;

  if (level_ <= 0) return;

  Index index = thisIndex;
  int icx,icy,icz;
  index.child(level_,&icx,&icy,&icz);
  index.set_level(level_ - 1);
  index.clean();

  // <duplicated code: refactor me!>
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();
  Restrict * restrict = simulation->problem()->restrict();

  FieldFace field_face (field_block,field_descr);

  field_face.set_restrict(restrict,icx,icy,icz);
  field_face.set_face(0,0,0); // 0-face is full block

  int narray; 
  char * array;
  field_face.load(&narray, &array);
  // </duplicated code>

  thisProxy[index].p_child_can_coarsen(icx,icy,icz,narray,array);

}

//----------------------------------------------------------------------

void CommBlock::p_child_can_coarsen(int icx,int icy, int icz,
				    int n, char * array)
{
  if (! can_coarsen()) return;

  // <duplicated code: refactor me!>
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  // allocate child block if this is the first
  if (!child_block_) {
    int nx,ny,nz;
    block_->field_block()->size(&nx,&ny,&nz);
    int num_field_blocks = block_->num_field_blocks();
    double xm,ym,zm;
    block_->lower(&xm,&ym,&zm);
    double xp,yp,zp;
    block_->upper(&xp,&yp,&zp);
    child_block_ = new Block  
      (nx, ny, nz, 
       num_field_blocks,
       xm, xp, ym, 
       yp, zm, zp);

    child_block_->allocate(field_descr);
  }

  FieldBlock * child_field_block = child_block_->field_block();
  FieldFace child_field_face(child_field_block, field_descr);
  Restrict *   restrict = simulation->problem()->restrict();
  child_field_face.set_restrict(restrict,icx,icy,icz);
  child_field_face.set_face(0,0,0);
   
  child_field_face.store (n, array);
  // </duplicated code>

  int rank = simulation->dimension();
  if (++count_coarsen_ >= NC(rank)) {
    delete block_;
    block_ = child_block_;
    child_block_ = 0;
    for (size_t i=0; i<children_.size(); i++) {
      if (children_[i] != thisIndex) {
	// Restricted child data is sent to parent when child destroyed
	thisProxy[children_[i]].ckDestroy();
      }
    }
  }
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_child (int n, char * buffer, 
				int icx, int icy, int icz)
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();
  Restrict * restrict = simulation->problem()->restrict();

  FieldFace field_face(field_block, field_descr);

  // Full block data
  field_face.set_face(0,0,0);

  field_face.set_restrict(restrict,icx,icy,icz);
   
  field_face.store (n, buffer);

  // Will need counter if not using QD
  // std::string refresh_type = simulation->config()->field_refresh_type;
  
  // if (refresh_type == "counter") {
  //   if (loop_refresh_.done()) {
  //     q_refresh_end();
  //   }
  // }
}

//----------------------------------------------------------------------

void CommBlock::q_adapt_stop()
{
  thisProxy.doneInserting();

  if (thisIndex.is_root()) {
    thisProxy.p_adapt_start();
  }
}

//----------------------------------------------------------------------

void CommBlock::q_adapt_end()
{
  TRACE("ADAPT CommBlock::q_adapt_end()");

  Performance * performance = simulation()->performance();
  if (performance->is_region_active(perf_adapt))
    performance->stop_region(perf_adapt);

  thisProxy.doneInserting();
  if (thisIndex.is_root()) {

    thisProxy.p_refresh_begin();

    // Check parameters now that all should have been accessed
    // (in particular Initial:foo:value, etc.)
    if (cycle() == simulation()->config()->initial_cycle) {
      simulation()->parameters()->check();
    }
  };



  TRACE ("END   PHASE ADAPT");
  
}

//----------------------------------------------------------------------

void CommBlock::set_child(Index index)
{ children_.push_back(index); }

//----------------------------------------------------------------------
  
#endif /* CONFIG_USE_CHARM */


