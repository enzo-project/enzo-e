// See LICENSE_CELLO file for license and copyright information

/// @file     charm_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions

#ifdef CONFIG_USE_CHARM
//#define DEBUG_ADAPT
#ifdef DEBUG_ADAPT

char buffer [80];

#define SET_FACE_LEVEL(INDEX,IF3,LEVEL,RECURSE,TYPE)		\
  sprintf (buffer,"set_face_level %d (%d %d %d  = %d) [%d]",	\
	   __LINE__,IF3[0],IF3[1],IF3[2],LEVEL,TYPE);		\
  INDEX.print(buffer,-1,2);					\
  thisProxy[INDEX].p_set_face_level (IF3,LEVEL,RECURSE,TYPE);
#else /* DEBUG_ADAPT */
#define SET_FACE_LEVEL(INDEX,IF3,LEVEL,RECURSE,TYPE)		\
  thisProxy[INDEX].p_set_face_level (IF3,LEVEL,RECURSE,TYPE);
#endif /* DEBUG_ADAPT */

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

const char * adapt_name[] = {
  "adapt_unknown", "adapt_same", "adapt_refine", "adapt_coarsen"
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
    adapt_step_ = initial_max_level;
  } else  {
    adapt_step_ = (mesh_max_level > 0) ? 1 : 0;
  }

  if (adapt_step_ > 0) {
    p_adapt_start();
  } else {
    q_adapt_end();
  }
}

//----------------------------------------------------------------------

void CommBlock::p_adapt_start()
{
  count_coarsen_ = 0;

  if (adapt_step_-- > 0) {

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
  TRACE("ADAPT CommBlock::q_adapt_next()");
  if (thisIndex.is_root()) {
    thisArray->doneInserting();
  }

  for (size_t i=0; i<child_face_level_.size(); i++) {
    child_face_level_[i] = face_level_unknown;
  }

  if (adapt_ == adapt_coarsen) coarsen();

  CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_stop(), 
			thisProxy[thisIndex]));
}

//----------------------------------------------------------------------

int CommBlock::determine_adapt()
{
  if (! is_leaf()) return adapt_same;

  FieldDescr * field_descr = simulation()->field_descr();

  int index_refine = 0;
  int adapt = adapt_unknown;

  Problem * problem = simulation()->problem();
  Refine * refine;

  while ((refine = problem->refine(index_refine++))) {

    int adapt_step_ = refine->apply(this, field_descr);

    adapt = reduce_adapt_(adapt,adapt_step_);

  }

  return adapt;
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

bool CommBlock::can_refine() const
{
  int max_level = simulation()->config()->mesh_max_level;
  return (is_leaf() && level_ < max_level);
}

//----------------------------------------------------------------------

void CommBlock::refine()
{
  if (! can_refine()) return;

  adapt_ = adapt_unknown;

  int rank = simulation()->dimension();
  
  int nc = NC(rank);

  int nx,ny,nz;
  block()->field_block()->size(&nx,&ny,&nz);

  int na3[3];
  size_forest (&na3[0],&na3[1],&na3[2]);

  int initial_cycle = simulation()->config()->initial_cycle;

  bool initial = (initial_cycle == cycle());

  for (int ic=0; ic<nc; ic++) {

    int ic3[3];
    ic3[0] = (ic & 1) >> 0;
    ic3[1] = (ic & 2) >> 1;
    ic3[2] = (ic & 4) >> 2;

    Index index_child = index_.index_child(ic3);

    if ( ! is_child(index_child) ) {

      int narray = 0;  
      char * array = 0;
      int iface[3] = {0,0,0};
      bool lghost[3] = {true,true,true};
      FieldFace * field_face = 
	load_face_ (&narray,&array, iface,ic3,lghost, op_array_prolong);

      const Factory * factory = simulation()->factory();

      int face_level[27];
      
      int if3[3];
      for (if3[0]=-1; if3[0]<=1; if3[0]++) {
	for (if3[1]=-1; if3[1]<=1; if3[1]++) {
	  for (if3[2]=-1; if3[2]<=1; if3[2]++) {
	    int ip3[3];
	    parent_face_(ip3,if3,ic3);
	    face_level[IF3(if3)] = face_level_[IF3(ip3)];
	  }
	}
      }
      
      int num_field_blocks = 1;
      bool testing = false;

      factory->create_block 
	(&thisProxy, index_child,
	 nx,ny,nz,
	 num_field_blocks,
	 adapt_step_,
	 initial,
	 cycle_,time_,dt_,
	 narray, array, op_array_prolong,
	 27,face_level,
	 testing);

      children_.push_back(index_child);

      delete field_face;

      refine_face_level_update_(index_child);

    }
  }
}

//----------------------------------------------------------------------

void CommBlock::refine_face_level_update_( Index index_child )
{
  Simulation * simulation = this->simulation();

  const int       rank         = simulation->dimension();
  const Config *  config       = simulation->config();
  const Problem * problem      = simulation->problem();
  const int       rank_refresh = config->field_refresh_rank;
  const bool      balance      = config->mesh_balance;
  const bool      periodic     = problem->boundary()->is_periodic();

  int na3[3];
  size_forest (&na3[0],&na3[1],&na3[2]);

  // Loop over all neighbors involved with communication

  int if3[3];
  ItFace it_face(rank,rank_refresh);
  while (it_face.next(if3)) {

    // MY NEIGHBOR

    Index index_neighbor = index_.index_neighbor
      (if3[0],if3[1],if3[2],na3,periodic);

    int jf3[3] = {-if3[0],-if3[1],-if3[2]};

    SET_FACE_LEVEL(index_neighbor,jf3,level_+1,false,adapt_refine);

    Index index_child_neighbor = index_child.index_neighbor
      (if3[0],if3[1],if3[2],na3,periodic);

    bool is_sibling = 
      index_child.index_parent() == index_child_neighbor.index_parent();

    if (is_sibling) {

      // CHILD SIBLING

      SET_FACE_LEVEL(index_child,if3,level_+1,false,adapt_refine);

      SET_FACE_LEVEL(index_child_neighbor,jf3,level_+1,false,adapt_refine);

    } else {

      int icc3[3];
      index_child.child(level_+1,icc3+0,icc3+1,icc3+2);
      int ip3[3];
      parent_face_(ip3,if3,icc3);

      int face_level = face_level_[IF3(ip3)];

      if (face_level == level_ - 1) {

	// CHILD GREAT UNCLE

	if (balance) {

	  // not level_-1 since balancing will force great uncle to be refined
	  SET_FACE_LEVEL(index_child,if3,level_,false,adapt_refine);

	  int ic3[3];
	  index_.child(level_,ic3+0,ic3+1,ic3+2);
	  int jc3[3];
	  facing_child_(jc3,ic3,if3);

	  Index index_uncle = index_neighbor.index_parent();

	  thisProxy[index_uncle].p_refine();

	}

      } else if (face_level == level_) {

	// CHILD UNCLE

	SET_FACE_LEVEL(index_child,if3,level_,false,adapt_refine);

	SET_FACE_LEVEL(index_neighbor,jf3,level_+1,true,adapt_refine);

      } else if (face_level == level_ + 1) {

	// CHILD COUSIN

	SET_FACE_LEVEL(index_child,if3,level_+1,false,adapt_refine);

	SET_FACE_LEVEL(index_child_neighbor,jf3,level_+1,false,adapt_refine);

      } else if (face_level == level_ + 2) {

	// CHILD NIBLINGS

	SET_FACE_LEVEL(index_child,if3,level_+2,false,adapt_refine);

	int ic3m[3],ic3p[3],ic3[3];
	loop_limits_nibling_(ic3m,ic3p, if3);

	for (ic3[0]=ic3m[0]; ic3[0]<=ic3p[0]; ic3[0]++) {
	  for (ic3[1]=ic3m[1]; ic3[1]<=ic3p[1]; ic3[1]++) {
	    for (ic3[2]=ic3m[2]; ic3[2]<=ic3p[2]; ic3[2]++) {

	      Index index_nibling = 
		index_child_neighbor.index_child(ic3[0],ic3[1],ic3[2]);

	      SET_FACE_LEVEL(index_nibling,jf3,level_+1,false,adapt_refine);
	    }
	  }
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void CommBlock::set_face_level (int if3[3], int level, int recurse, int type)
{ 

  int index_face = IF3(if3);

  if (face_level_[index_face] == face_level_unknown) {
    face_level_[index_face] = level;
  } else if (type == adapt_refine) {
    face_level_[index_face] = std::max(face_level_[index_face],level); 
  } else if (type == adapt_coarsen) {
    face_level_[index_face] = std::min(face_level_[index_face],level); 
  }

  if (recurse) {

    if (type == adapt_refine && ! is_leaf()) {

      // update children that share the same face

      for (size_t ic=0; ic<children_.size(); ic++) {

	Index index_child = children_[ic];

	if (index_child != index_) {

	  int ic3[3];
	  index_child.child(level_+1,ic3,ic3+1,ic3+2);

	  bool face_adjacent = child_is_on_face_(ic3,if3);

	  if (face_adjacent) {

	    SET_FACE_LEVEL(index_child,if3,level,false,adapt_refine);

	    int ifc3m[3], ifc3p[3], ifc3[3];
	    loop_limits_faces_(ifc3m,ifc3p,if3,ic3);

	    for (ifc3[0]=ifc3m[0]; ifc3[0]<=ifc3p[0]; ifc3[0]++) {
	      for (ifc3[1]=ifc3m[1]; ifc3[1]<=ifc3p[1]; ifc3[1]++) {
		for (ifc3[2]=ifc3m[2]; ifc3[2]<=ifc3p[2]; ifc3[2]++) {
		  SET_FACE_LEVEL(index_child,ifc3,level,false,adapt_refine);
		}
	      }
	    }
	  }
	}
      }
    } else if (type == adapt_coarsen && level_ > 0) {
      Index index_parent = index_.index_parent();
      int ic3[3];
      index_.child(level_,ic3+0,ic3+1,ic3+2);
      int ip3[3];
      parent_face_(ip3,if3,ic3);
      SET_FACE_LEVEL(index_parent,ip3,level,false,adapt_coarsen);
    }

  } 
}

//----------------------------------------------------------------------

bool CommBlock::child_is_on_face_(int ic3[3], int if3[3]) const
{
  bool face_adjacent = true;
  if (if3[0] == -1 && ic3[0] != 0) face_adjacent = false;
  if (if3[0] == +1 && ic3[0] != 1) face_adjacent = false;
  if (if3[1] == -1 && ic3[1] != 0) face_adjacent = false;
  if (if3[1] == +1 && ic3[1] != 1) face_adjacent = false;
  if (if3[2] == -1 && ic3[2] != 0) face_adjacent = false;
  if (if3[2] == +1 && ic3[2] != 1) face_adjacent = false;
  return face_adjacent;
}

//----------------------------------------------------------------------

void CommBlock::loop_limits_faces_ 
(int ic3m[3], int ic3p[3], int if3[3], int ic3[3]) const
{
  ic3m[0] = (if3[0] != 0) ? if3[0] :  -ic3[0];
  ic3p[0] = (if3[0] != 0) ? if3[0] : 1-ic3[0];
  ic3m[1] = (if3[1] != 0) ? if3[1] :  -ic3[1];
  ic3p[1] = (if3[1] != 0) ? if3[1] : 1-ic3[1];
  ic3m[2] = (if3[2] != 0) ? if3[2] :  -ic3[2];
  ic3p[2] = (if3[2] != 0) ? if3[2] : 1-ic3[2];
	
  Simulation * simulation = this->simulation();
  const int rank = simulation->dimension();

  if (rank < 2) { ic3m[1]=0; ic3p[1]=0; }
  if (rank < 3) { ic3m[2]=0; ic3p[2]=0; }
}

//----------------------------------------------------------------------

void CommBlock::parent_face_(int ip3[3],int if3[3], int ic3[3]) const
{
  ip3[0] = if3[0];
  ip3[1] = if3[1];
  ip3[2] = if3[2];
  if (if3[0] == +1 && ic3[0] == 0) ip3[0] = 0;
  if (if3[0] == -1 && ic3[0] == 1) ip3[0] = 0;
  if (if3[1] == +1 && ic3[1] == 0) ip3[1] = 0;
  if (if3[1] == -1 && ic3[1] == 1) ip3[1] = 0;
  if (if3[2] == +1 && ic3[2] == 0) ip3[2] = 0;
  if (if3[2] == -1 && ic3[2] == 1) ip3[2] = 0;
}

//----------------------------------------------------------------------

bool CommBlock::can_coarsen() const
{ 
  if (level_ <= 0) return false;

  if (! is_leaf() ) return false;

  int refresh_rank = simulation()->config()->field_refresh_rank;
  int rank = simulation()->dimension();
  int if3[3];

  ItFace it_face(rank,refresh_rank);

  while (it_face.next(if3)) {
    int i=IF3(if3);
    if (face_level_[i] > level_) return false;
  }

  return true;
}

//----------------------------------------------------------------------

void CommBlock::coarsen()
{
  if (! can_coarsen()) return;
  
  adapt_ = adapt_unknown;

  Index index_parent = thisIndex;
  int ichild[3];
  index_parent.child(level_,ichild,ichild+1,ichild+2);
  index_parent.set_level(level_ - 1);
  index_parent.clean();

  int narray; 
  char * array;
  int iface[3] = {0,0,0};
  bool lghost[3] = {true,true,true};
  FieldFace * field_face = 
    load_face_(&narray,&array, iface, ichild, lghost, op_array_restrict);

  int nf = face_level_.size();
  int face_level[nf];
  for (int i=0; i<nf; i++) face_level[i] = face_level_[i];

  thisProxy[index_parent].p_child_can_coarsen
    (ichild, narray,array, nf,face_level);

  delete field_face;

}

//----------------------------------------------------------------------

FieldFace * CommBlock::load_face_
(
 int * narray, char ** array, 
 int iface[3], int ichild[3], bool lghost[3],
 int op_array_type
)
{
  FieldFace * field_face = create_face_ (iface,ichild,lghost,
					 op_array_type);

  field_face->load(narray, array);
  return field_face;
}

//----------------------------------------------------------------------

FieldFace * CommBlock::store_face_
(
 int narray, char * array, 
 int iface[3], int ichild[3], bool lghost[3],
 int op_array_type
)
{
  FieldFace * field_face = create_face_ (iface,ichild,lghost,
					 op_array_type);

  field_face->store(narray, array);
  delete field_face;
  return NULL;
}

//----------------------------------------------------------------------

FieldFace * CommBlock::create_face_
(int iface[3], int ichild[3], bool lghost[3],
 int op_array_type
)
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();

  FieldFace * field_face = new FieldFace (field_block,field_descr);

  if (op_array_type == op_array_restrict) {
    Restrict * restrict = simulation->problem()->restrict();
    field_face->set_restrict(restrict,ichild[0],ichild[1],ichild[2]);
  } else if (op_array_type == op_array_prolong) {
    Prolong * prolong = simulation->problem()->prolong();
    field_face->set_prolong(prolong,ichild[0],ichild[1],ichild[2]);
  }

  field_face->set_face(iface[0],iface[1],iface[2]);
  field_face->set_ghost(lghost[0],lghost[1],lghost[2]);
  return field_face;
}

//----------------------------------------------------------------------

void CommBlock::p_child_can_coarsen(int ichild[3],
				    int na, char * array,
				    int nf, int * child_face_level)
{
#ifdef COARSEN
#else /* COARSEN */
  return;
#endif /* COARSEN */

  // allocate child block if this is the first
  if (!child_block_) {

    FieldDescr * field_descr = simulation()->field_descr();
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

  int iface[3] = {0,0,0};
  bool lghost[3] = {true,true,true};
  store_face_(na,array, iface, ichild, lghost, op_array_restrict);

  int rank = simulation()->dimension();
  int refresh_rank = simulation()->config()->field_refresh_rank;

  int if3[3];
  ItFace it_face (rank,refresh_rank);

  // update child_is_on_face_[] with subset of child's
  while (it_face.next(if3)) {
    if (child_is_on_face_(ichild,if3)) {
      int index_face = IF3(if3);
      child_face_level_[index_face] = 
	std::min (child_face_level_[index_face],
		  child_face_level[index_face]);
    }
  }

  if (++count_coarsen_ >= NC(rank)) {

    delete block_;
    block_ = child_block_;
    child_block_ = 0;

    // need to update application so it knows where fields arrays are
    initialize ();

    face_level_ = child_face_level_;

    for (size_t i=0; i<children_.size(); i++) {
      Index index_child = children_[i];
      thisProxy[index_child].p_parent_coarsened();
      coarsen_face_level_update_(index_child);
    }
    children_.clear();

  }
}

//----------------------------------------------------------------------

void CommBlock::p_parent_coarsened()
{
  coarsened_ = true;
}

//----------------------------------------------------------------------

void CommBlock::coarsen_face_level_update_( Index index_child )
{
  Simulation * simulation = this->simulation();
  const Config *  config  = simulation->config();
  const Problem * problem = simulation->problem();

  const int  rank         = simulation->dimension();
  const int  rank_refresh = config->field_refresh_rank;
  const bool periodic     = problem->boundary()->is_periodic();

  int na3[3];
  size_forest (&na3[0],&na3[1],&na3[2]);

  // Loop over all neighbors involved with communication

  ItFace it_face(rank,rank_refresh);
  int if3[3];
  while (it_face.next(if3)) {

    int jf3[3] = {-if3[0],-if3[1],-if3[2]};

    Index index_neighbor = index_.index_neighbor
      (if3[0],if3[1],if3[2],na3,periodic);

    Index index_child_neighbor = index_child.index_neighbor
      (if3[0],if3[1],if3[2],na3,periodic);

    bool is_sibling = 
      index_child.index_parent() == index_child_neighbor.index_parent();

    if (is_sibling) {

    } else {

      int icc3[3];
      index_child.child(level_+1,icc3+0,icc3+1,icc3+2);
      int ip3[3];
      parent_face_(ip3,if3,icc3);

      int face_level = face_level_[IF3(ip3)];

      if (face_level == level_ - 1) {
	
	Index index_uncle = index_neighbor.index_parent();

	SET_FACE_LEVEL(index_uncle,jf3,level_,false,adapt_coarsen);

      } else if (face_level == level_) {

 	SET_FACE_LEVEL(index_neighbor,jf3,level_,false,adapt_coarsen);
	
      } else if (face_level == level_ + 1) {

	SET_FACE_LEVEL(index_child_neighbor,jf3,level_,true,adapt_coarsen);

	int ic3m[3],ic3p[3],ic3[3];
	loop_limits_nibling_(ic3m,ic3p, if3);

	for (ic3[0]=ic3m[0]; ic3[0]<=ic3p[0]; ic3[0]++) {
	  for (ic3[1]=ic3m[1]; ic3[1]<=ic3p[1]; ic3[1]++) {
	    for (ic3[2]=ic3m[2]; ic3[2]<=ic3p[2]; ic3[2]++) {

	      Index index_nibling = 
		index_child_neighbor.index_child(ic3[0],ic3[1],ic3[2]);

	      SET_FACE_LEVEL(index_nibling,jf3,level_+1,false,adapt_coarsen);
	    }
	  }
	}
      } else {
	index_.print("index_");
	index_child.print("index_child");
	ERROR6("CommBlock::coarsen_face_level_update_()",
	       "Unhandled level jump from level_[%d] to face_level(%d %d %d)[%d] = %d",
	       level_,ip3[0],ip3[1],ip3[2],IF3(ip3),face_level);
      }
    }
  }
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_child (int n, char * buffer, int ichild[3])
{
  int iface[3] = {0,0,0};
  bool lghost[3] = {true,true,true};
  store_face_(n,buffer, iface, ichild, lghost, op_array_restrict);
}

//----------------------------------------------------------------------

void CommBlock::q_adapt_stop()
{
  //  TRACE_CHARM("doneInserting");
  thisProxy.doneInserting();
  TRACE("ADAPT CommBlock::q_adapt_stop()");
  if (thisIndex.is_root()) {
    thisArray->doneInserting();
    thisProxy.p_adapt_start();
  }
}

//----------------------------------------------------------------------

void CommBlock::q_adapt_end()
{
  // CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_exit(), 
  // 			thisProxy[thisIndex]));

  //  thisProxy.doneInserting();
  TRACE("ADAPT CommBlock::q_adapt_end()");

  debug_faces_("child",&face_level_[0]);

  TRACE ("END   PHASE ADAPT\n");

  next_phase_ = phase_output;

  if (coarsened_) {
    TRACE_CHARM("ckDestroy()");
    thisProxy[thisIndex].ckDestroy();
    thisProxy.doneInserting();
    return;
  } else {
    //    TRACE_CHARM("doneInserting");
    thisProxy.doneInserting();
  }

  CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_exit(), 
			thisProxy[thisIndex]));
}

//----------------------------------------------------------------------

void CommBlock::q_adapt_exit()
{

  Performance * performance = simulation()->performance();
  if (performance->is_region_active(perf_adapt)) {
    performance->stop_region(perf_adapt);
  }

  if (thisIndex.is_root()) {


    thisProxy.p_refresh_begin();

    // Check parameters now that all should have been accessed
    // (in particular Initial:foo:value, etc.)
    if (cycle() == simulation()->config()->initial_cycle) {
      simulation()->parameters()->check();
    }
  }
}

//----------------------------------------------------------------------
  
#endif /* CONFIG_USE_CHARM */


