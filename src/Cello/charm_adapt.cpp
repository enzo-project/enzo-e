// See LICENSE_CELLO file for license and copyright information

/// @file     charm_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions

#ifdef CONFIG_USE_CHARM

/* #define DEBUG_ADAPT */

#ifdef DEBUG_ADAPT

char buffer [80];

#define SET_FACE_LEVEL(INDEX,IF3,LEVEL,RECURSE,LINE)		\
  sprintf (buffer,"set_face_level(%d %d %d  = %d) [%d]",	\
	   IF3[0],IF3[1],IF3[2],LEVEL,__LINE__);		\
  INDEX.print(buffer,-1,2);					\
  thisProxy[INDEX].p_set_face_level (IF3,LEVEL,RECURSE,LINE);
#else /* DEBUG_ADAPT */
#define SET_FACE_LEVEL(INDEX,IF3,LEVEL,RECURSE,LINE)		\
  thisProxy[INDEX].p_set_face_level (IF3,LEVEL,RECURSE,LINE);
#endif /* DEBUG_ADAPT */

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
  TRACE("ADAPT CommBlock::q_adapt_next()");
  thisProxy.doneInserting();

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

  bool initial = (initial_cycle == cycle());

  for (int ic=0; ic<nc; ic++) {

    int ic3[3];
    ic3[0] = (ic & 1) >> 0;
    ic3[1] = (ic & 2) >> 1;
    ic3[2] = (ic & 4) >> 2;

    Index index_child = index_.index_child(ic3);

    if ( ! is_child(index_child) ) {

      // create child

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
      
#ifdef FULL_GHOST
      field_face.set_ghost(true,true,true);
#endif

      field_face.load(&narray, &array);

      //       if (cycle_ > 0) index_.print("refine",-1,2,true);
	

      // </duplicated code>

      const Factory * factory = simulation->factory();

      int face_level[27];
      
      int if3[3],if3m[3],if3p[3];
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

      //       if (cycle_ > 0) printf ("load %d\n",narray);

      factory->create_block 
	(&thisProxy, index_child,
	 nx,ny,nz,
	 num_field_blocks,
	 count_adapt_,
	 initial,
	 cycle_,time_,dt_,
	 narray, array, op_array,
	 27,face_level,
	 testing);

      set_child(index_child);

      // update child and neighbor face_level[]

      face_level_update_new_(index_child);

    }
  }
}

//----------------------------------------------------------------------

void CommBlock::face_level_update_new_( Index index_child )
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

  int if3m[3], if3p[3], if3[3];
  loop_limits_refresh_(if3m+0,if3m+1,if3m+2,
		       if3p+0,if3p+1,if3p+2);

  if3m[0] = (rank >= 1) ? -1 : 0;
  if3m[1] = (rank >= 2) ? -1 : 0;
  if3m[2] = (rank >= 3) ? -1 : 0;
  if3p[0] = (rank >= 1) ? +1 : 0;
  if3p[1] = (rank >= 2) ? +1 : 0;
  if3p[2] = (rank >= 3) ? +1 : 0;

  // Loop over all neighbors involved with communication

  for (if3[0]=if3m[0]; if3[0]<=if3p[0]; if3[0]++) {
    for (if3[1]=if3m[1]; if3[1]<=if3p[1]; if3[1]++) {
      for (if3[2]=if3m[2]; if3[2]<=if3p[2]; if3[2]++) {

	int rank_face = 
	  rank - (abs(if3[0]) + abs(if3[1]) + abs(if3[2]));

	if (rank_refresh <= rank_face && rank_face < rank) {

	  // MY NEIGHBOR

	  Index index_neighbor = index_.index_neighbor
	    (if3[0],if3[1],if3[2],na3,periodic);

	  int jf3[3] = {-if3[0],-if3[1],-if3[2]};

	  SET_FACE_LEVEL(index_neighbor,jf3,level_+1,false,__LINE__);

	  // CHILD NEIGHBOR

	  Index index_child_neighbor = index_child.index_neighbor
	    (if3[0],if3[1],if3[2],na3,periodic);

	  bool is_sibling = 
	    index_child.index_parent() == index_child_neighbor.index_parent();

	  if (is_sibling) {

	    // SIBLING

	    SET_FACE_LEVEL(index_child,if3,level_+1,false,__LINE__);

	    SET_FACE_LEVEL(index_child_neighbor,jf3,level_+1,false,__LINE__);

	  } else {

	    int icc3[3];
	    index_child.child(level_+1,icc3+0,icc3+1,icc3+2);
	    int ip3[3];
	    parent_face_(ip3,if3,icc3);

	    int face_level = face_level_[IF3(ip3)];

	    if (face_level == level_ - 1) {

	      // GREAT UNCLE

	      if (balance) {
		Index index_uncle = index_neighbor.index_parent();

		SET_FACE_LEVEL(index_child,if3,level_,false,__LINE__);

		// facing child
		int ic3[3];
		index_.child(level_,ic3+0,ic3+1,ic3+2);
		int jc3[3] = {ic3[0],ic3[1],ic3[2]};
		if (if3[0] == -1) jc3[0] = 1;
		if (if3[1] == -1) jc3[1] = 1;
		if (if3[2] == -1) jc3[2] = 1;
		if (if3[0] == +1) jc3[0] = 0;
		if (if3[1] == +1) jc3[1] = 0;
		if (if3[2] == +1) jc3[2] = 0;

		thisProxy[index_uncle].p_balance (jc3,jf3,level_+1);

	      }

	    } else if (face_level == level_) {

	      // UNCLE

	      SET_FACE_LEVEL(index_child,if3,level_,false,__LINE__);

	      SET_FACE_LEVEL(index_neighbor,jf3,level_+1,true,__LINE__);

	    } else if (face_level == level_ + 1) {

	      // COUSIN

	      SET_FACE_LEVEL(index_child,if3,level_+1,false,__LINE__);

	      SET_FACE_LEVEL(index_child_neighbor,jf3,level_+1,false,__LINE__);

	    } else if (face_level == level_ + 2) {

	      // NIBLINGS

	      SET_FACE_LEVEL(index_child,if3,level_+2,false,__LINE__);

	      int ic3m[3],ic3p[3],ic3[3];
	      loop_limits_nibling_(ic3m,ic3m+1,ic3m+2,
				   ic3p,ic3p+1,ic3p+2,
				   if3[0],if3[1],if3[2]);

	      for (ic3[0]=ic3m[0]; ic3[0]<=ic3p[0]; ic3[0]++) {
		for (ic3[1]=ic3m[1]; ic3[1]<=ic3p[1]; ic3[1]++) {
		  for (ic3[2]=ic3m[2]; ic3[2]<=ic3p[2]; ic3[2]++) {

		    Index index_nibling = 
		      index_child_neighbor.index_child(ic3[0],ic3[1],ic3[2]);

		    SET_FACE_LEVEL(index_nibling,jf3,level_+1,false,__LINE__);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


void CommBlock::set_face_level (int if3[3], int level, int recurse, int line)
{ 
  // recurse = false;
#ifdef DEBUG_ADAPT
  char buffer[80];
  sprintf (buffer,"SET_FACE_LEVEL(%d %d %d = %d [%d]",if3[0],if3[1],if3[2],level,line);
  index_.print(buffer,-1,2);
#endif

  face_level_[IF3(if3)] = std::max(face_level_[IF3(if3)],level); 
  //  face_level_[IF3(if3)] = level;

  if (recurse && ! is_leaf()) {

    // Simulation * simulation = this->simulation();
    // const int rank         = simulation->dimension();
    // const int rank_refresh = config->field_refresh_rank;
    // const int rank_face = abs(if3[0]) + abs(if3[1]) + abs(if3[2]);

    // update children that share the same face

    for (size_t ic=0; ic<children_.size(); ic++) {

      Index index_child = children_[ic];

      if (index_child != index_) {

	int ic3[3];
	index_child.child(level_+1,ic3,ic3+1,ic3+2);

	bool face_adjacent = child_is_on_face_(if3,ic3);

	if (face_adjacent) {

#ifdef CONFIG_USE_CHARM
	  SET_FACE_LEVEL(index_child,if3,level,true,__LINE__);
#endif

	  int ifc3m[3], ifc3p[3], ifc3[3];
	  loop_limits_faces_(ifc3m,ifc3p,if3,ic3);

	  for (ifc3[0]=ifc3m[0]; ifc3[0]<=ifc3p[0]; ifc3[0]++) {
	    for (ifc3[1]=ifc3m[1]; ifc3[1]<=ifc3p[1]; ifc3[1]++) {
	      for (ifc3[2]=ifc3m[2]; ifc3[2]<=ifc3p[2]; ifc3[2]++) {
#ifdef CONFIG_USE_CHARM
		SET_FACE_LEVEL(index_child,ifc3,level,true,__LINE__);
#endif
	      }
	    }
	  }
	}
      }
    }
  } 
}

//----------------------------------------------------------------------

bool CommBlock::child_is_on_face_(int if3[3], int ic3[3]) const
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
  ic3m[0] = (if3[0] != 0) ? if3[0] : -ic3[0];
  ic3p[0] = (if3[0] != 0) ? if3[0] : 1-ic3[0];
  ic3m[1] = (if3[1] != 0) ? if3[1] : -ic3[1];
  ic3p[1] = (if3[1] != 0) ? if3[1] : 1-ic3[1];
  ic3m[2] = (if3[2] != 0) ? if3[2] : -ic3[2];
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

void CommBlock::p_balance(int ic3[3], int if3[3], int level)
{
  refine();
  //  Index index_child = index_.index_child(ic3[0],ic3[1],ic3[2]);
  //  SET_FACE_LEVEL(index_child,if3,level,false,__LINE__);
}

//----------------------------------------------------------------------

bool CommBlock::can_coarsen() const
{ 
  if (level_ <= 0) return false;

  int refresh_rank = simulation()->config()->field_refresh_rank;
  int rank = simulation()->dimension();
  int if3m[3],if3p[3],if3[3];

  // loop_limits_refresh_(if3m+0,if3m+1,if3m+2,
  // 		       if3p+0,if3p+1,if3p+2);
  if3m[0] = (rank >= 1) ? -1 : 0;
  if3m[1] = (rank >= 2) ? -1 : 0;
  if3m[2] = (rank >= 3) ? -1 : 0;
  if3p[0] = (rank >= 1) ? +1 : 0;
  if3p[1] = (rank >= 2) ? +1 : 0;
  if3p[2] = (rank >= 3) ? +1 : 0;

  for (if3[0]=if3m[0]; if3[0]<=if3p[0]; if3[0]++) {
    for (if3[1]=if3m[1]; if3[1]<=if3p[1]; if3[1]++) {
      for (if3[2]=if3m[2]; if3[2]<=if3p[2]; if3[2]++) {
	int face_rank = rank - (abs(if3[0]) + abs(if3[1]) + abs(if3[2]));
	if (refresh_rank <= face_rank && face_rank < rank) {
	  int i=IF3(if3);
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
  if (! can_coarsen()) return;
  
  adapt_ = adapt_unknown;

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
#ifdef FULL_GHOST
  field_face.set_ghost(true,true,true);
#endif

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
#ifdef FULL_GHOST
  child_field_face.set_ghost(true,true,true);
#endif
   
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
#ifdef FULL_GHOST
  field_face.set_ghost(true,true,true);
#endif

  field_face.set_restrict(restrict,icx,icy,icz);
   
  field_face.store (n, buffer);

}

//----------------------------------------------------------------------

void CommBlock::q_adapt_stop()
{
  TRACE("ADAPT CommBlock::q_adapt_stop()");
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

  debug_faces_("child",&face_level_[0]);

  TRACE ("END   PHASE ADAPT\n");

  next_phase_ = phase_output;

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

void CommBlock::set_child(Index index)
{ children_.push_back(index); }

//----------------------------------------------------------------------
  
#endif /* CONFIG_USE_CHARM */


