// See LICENSE_CELLO file for license and copyright information

#ifdef TEMP_USE_COMM

/// @file     comm_CommBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-11-18
/// @brief    Communication class associated with Blocks

#include "comm.hpp"

#ifdef CONFIG_USE_CHARM

//----------------------------------------------------------------------

CommBlock::CommBlock
(
 int ibx, int iby, int ibz,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xpm, double ypm, double zpm, // Patch begin
 double xb, double yb, double zb,    // Block width
 CProxy_Patch proxy_patch,
 int patch_id,
 int patch_rank,
 int num_field_blocks
) throw ()
  :  count_refresh_face_(0)
{
  block_ = new Block 
    (ibx,iby,ibz, 
     nbx,nby,nbz,
     nx,ny,nz,
     xm,ym,zm, 
     xb,yb,zb, 
     patch_id,
     patch_rank,
     num_field_blocks,
     this);
}

//----------------------------------------------------------------------

CommBlock::CommBlock
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xpm, double ypm, double zpm, // Patch begin
 double xb, double yb, double zb,    // Block width
 CProxy_Patch proxy_patch,
 int patch_id,
 int patch_rank,
 int num_field_blocks) throw ()
{
  block = new Block
    (nbx,nby,nbz,
     nx,ny,nz,
     xm,ym,zm, 
     xb,yb,zb, 
     proxy_patch,
     patch_id,
     patch_rank,
     num_field_blocks,
     this);
}

//----------------------------------------------------------------------

CommBlock::~CommBlock() throw ()
{
}

//----------------------------------------------------------------------

CommBlock::CommBlock(const CommBlock & CommBlock) throw ()
/// @param     CommBlock  Object being copied
{
}

//----------------------------------------------------------------------

CommBlock & CommBlock::operator= (const CommBlock & CommBlock) throw ()
/// @param     CommBlock  Source object of the assignment
/// @return    The target assigned object
{
  return *this;
}

//----------------------------------------------------------------------

void CommBlock::p_initial()
{
  Simulation * simulation  = proxy_simulation.ckLocalBranch();
  block_->initial(simulation);

  proxy_patch_.s_initial();

}

//----------------------------------------------------------------------

void CommBlock::p_refresh() 
{
  block_->refresh(); 
}

//----------------------------------------------------------------------

void CommBlock::p_compute(int cycle, double time, double dt)
{
  block_->compute(); 
}

//----------------------------------------------------------------------

void CommBlock::p_output(CkReductionMsg * msg)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  block_->output(simulation);
  proxy_patch_.s_output(Simulation * simulation);

}

//----------------------------------------------------------------------

void CommBlock::p_read (int index_initial)
{
  INCOMPLETE("CommBlock::p_read");
}

//----------------------------------------------------------------------

void CommBlock::p_exchange (int n, char * buffer, int fx, int fy, int fz)
{

  TRACE ("CommBlock::p_exchange()");
  
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  // Call exchange() on the Block
  block_->exchange(simulation,int n, char * buffer, int fx, int fy, int fz);

  // When done exchanging, call prepare

  if (++count_refresh_face_ >= count) {
    count_refresh_face_ = 0;
    prepare();
  }
}

//----------------------------------------------------------------------

void CommBlock::p_write (int index_output)
{
  TRACE("OUTPUT CommBlock::p_write()");

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  // Call output() on the Block

  FieldDescr * field_descr = simulation->field_descr();
  Problem    * problem     = simulation->problem();
  Output     * output      = problem->output(index_output);

  output->write_block(block_,field_descr,0,0,0);

  proxy_patch_.s_write();
}

//======================================================================

#endif /* CONFIG_USE_CHARM */

#endif /* TEMP_USE_COMM */
