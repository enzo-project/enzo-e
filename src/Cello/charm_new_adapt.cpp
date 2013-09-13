// See LICENSE_CELLO file for license and copyright information

/// @file     charm_new_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions
///
/// create_mesh()
///
/// 1. Apply refinement criteria
/// 2. Notify neighbors of intent
/// 3. Update intent based on neighbor's (goto 2)
/// 4. After quiescence, insert() / delete() array elements
/// 5. After contribute(), doneInserting()

/// adapt_mesh()
///
/// 1. Apply refinement criteria
/// 2. Notify neighbors of intent
/// 3. Update intent based on neighbor's (goto 2)
/// 4. After quiescence, insert() / delete() array elements
///
/// 1. Apply refinement criteria
///
///    adapt_mesh() or create_mesh()
///
///    Call determine_adapt() to decide whether each existing leaf
///    block should COARSEN, REFINE or stay the SAME.  If REFINE,
///    create new child blocks, including upsampled data if refining.
///    If creating initial mesh, the new CommBlocks will apply initial
///    conditions, and recursively call create_mesh.  After
///    quiescence, proceed to next AMR step.
///
/// 2. Notify all neighbors of intent
///
///    notify_neighbors()
///
///    If REFINE or SAME, loop over all neighbors and call
///    index_neighbor.p_set_face_level() to notify them of intended
///    level.
///
/// 3. Update intent based on neighbor's (goto 2)
///
///    p_set_face_level() updates intended level based on neighbors
///    intended level.  If intended level changes, call
///    notify_neighbors() to notify all neighbors of updated intent
///
/// 4. After quiescence, insert() / delete() array elements
///
///    q_adapt_finalize()
///
///    After quiescence, all refinement decisions are known.  All
///    REFINE elements should already be created, so only COARSEN
///    needs to be called on elements that are to be coarsened.  Since
///    all REFINE elements are created, a call to doneInserting() is
///    performed on the chare array.
///
///----------------------------------------------------------------------

/* #define CELLO_TRACE */

#ifdef TEMP_NEW_ADAPT

//--------------------------------------------------
#ifdef DEBUG_ADAPT

char buffer [80];

#define SET_FACE_LEVEL(INDEX,IF3,LEVEL)		\
  sprintf (buffer,"p_get_neighbor_level %d (%d  = %d) [%d]",	\
	   __LINE__,IF3[0],IF3[1],IF3[2]);		\
  INDEX.print(buffer,-1,2);					\
  thisProxy[INDEX].p_get_neighbor_level (IF3,LEVEL);
#else /* DEBUG_ADAPT */
#define SET_FACE_LEVEL(INDEX,IF3,LEVEL)		\
  thisProxy[INDEX].p_get_neighbor_level (IF3,LEVEL);

#endif /* DEBUG_ADAPT */
//--------------------------------------------------

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

const char * adapt_name[] = {
  "adapt_unknown", "adapt_same", "adapt_refine", "adapt_coarsen"
};

//======================================================================

void CommBlock::create_mesh()
{
  TRACE("create_mesh()");
  int level_maximum = simulation()->config()->initial_max_level;

  level_desired_ = desired_level_(level_maximum);

  notify_neighbors(level_desired_);

  CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_end(), 
			thisProxy[thisIndex]));

}

//----------------------------------------------------------------------

void CommBlock::adapt_mesh()
{
  TRACE("adapt_mesh()");
  int level_maximum = simulation()->config()->mesh_max_level;

  level_desired_ = desired_level_(level_maximum);

  notify_neighbors(level_desired_);

  CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_end(), 
			thisProxy[thisIndex]));

}

//----------------------------------------------------------------------

int CommBlock::desired_level_(int level_maximum)
{
  int level = this->level();
  adapt_ = determine_adapt();
  if (adapt_ == adapt_coarsen && level > 0) 
    --level;
  if (adapt_ == adapt_refine  && level < level_maximum) 
      ++level;
  return level;
}

//----------------------------------------------------------------------

void CommBlock::q_adapt_end()
{
  TRACE0;
  next_phase_ = phase_output;
  if (thisIndex.is_root()) {

    thisProxy.p_refresh_begin();
  }
}

//----------------------------------------------------------------------

void CommBlock::notify_neighbors(int level)
{
  // Loop over all neighobrs and (if not coarsening) tell them desired refinement level
  if (adapt_ != adapt_coarsen) {

    Simulation * simulation = this->simulation();
    Problem * problem       = simulation->problem();
    const Config *  config  = simulation->config();
    const int  rank         = simulation->dimension();
    const int  rank_refresh = config->field_refresh_rank;
    const bool periodic     = problem->boundary()->is_periodic();

    int na3[3];
    size_forest (&na3[0],&na3[1],&na3[2]);

    ItFace it_face(rank,rank_refresh);
    int if3[3];
    while (it_face.next(if3)) {
      Index index_neighbor = index_.index_neighbor
	(if3[0],if3[1],if3[2],na3,periodic);
      SET_FACE_LEVEL(index_neighbor,if3,level);
    }
  }
}

//----------------------------------------------------------------------

void CommBlock::p_get_neighbor_level(int if3[3], int level)
{
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

void CommBlock::x_refresh_child (int n, char * buffer, int ic3[3])
{
  int  iface[3]  = {0,0,0};
  bool lghost[3] = {true,true,true};
  store_face_(n,buffer, iface, ic3, lghost, op_array_restrict);
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

#endif /* #ifdef TEMP_NEW_ADAPT */

