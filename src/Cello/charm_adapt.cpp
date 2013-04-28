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

#define IC(ix,iy,iz)  ((ix+2)%2) + 2*( ((iy+2)%2) + 2*((iz+2)%2) )
#define NC(array_rank) (array_rank == 1 ? 2 : (array_rank == 2 ? 4 : 8))
#define IN(axis,face)  ((face) + 2*(axis))
#define NN(array_rank) (array_rank*2)

const char * adapt_name[] = {
  "adapt_unknown",
  "adapt_same",
  "adapt_refine",
  "adapt_coarsen"
};

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
  int adapt = adapt_unknown;
  while (Refine * refine = simulation_charm->problem()->refine(i++)) {
    adapt = update_adapt_(adapt,refine->apply(field_block, field_descr));
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

void CommBlock::refine()
{
  TRACE("ADAPT CommBlock::refine()");

  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();
  const Factory * factory = simulation_charm->factory();
  int rank = simulation_charm->dimension();
  int nc = NC(rank);

  int ibx = index_[0];
  int iby = index_[1];
  int ibz = index_[2];
  int nbx = size_[0];
  int nby = size_[1];
  int nbz = size_[2];

  int nx,ny,nz;
  block()->field_block()->size(&nx,&ny,&nz);

  double xm[2],ym[2],zm[2];
  block()->lower(&xm[0],&ym[0],&zm[0]);
  double xp[2],yp[2],zp[2];
  block()->upper(&xm[1],&ym[1],&zm[1]);
  xm[1] = xp[0] = 0.5 * (xm[0] + xp[1]);
  ym[1] = yp[0] = 0.5 * (ym[0] + yp[1]);
  zm[1] = zp[0] = 0.5 * (zm[0] + zp[1]);

  int num_field_blocks = 1;
  bool testing = false;

  for (int ic=0; ic<nc; ic++) {

    int ix = (ic & 1) >> 0;
    int iy = (ic & 2) >> 1;
    int iz = (ic & 4) >> 2;

    TRACE5("ADAPT new child %d [%d %d %d]/ %d",ic,ix,iy,iz,nc);

    Index index = thisIndex;
    index.set_level(level_+1);
    index.set_tree (level_+1,ix,iy,iz);
    index.clean();

    factory->create_block 
      (thisProxy, &index,
       ibx,iby,ibz,
       nbx,nby,nbz,
       nx,ny,nz,
       level_+1,
       xm[ix],ym[iy],zm[iz],
       xp[ix],yp[iy],zp[iz],
       num_field_blocks,
       testing);

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
  TRACE1("ADAPT q_adapt(%d)",level_active_);
  ++level_active_;
  adapt();
}

//----------------------------------------------------------------------

void CommBlock::p_balance()
{
  TRACE("ADAPT CommBlock::p_balance()");
}

//----------------------------------------------------------------------


//======================================================================

#endif /* CONFIG_USE_CHARM */


