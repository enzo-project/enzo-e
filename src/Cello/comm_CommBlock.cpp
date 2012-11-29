// See LICENSE_CELLO file for license and copyright information

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

void CommBlock::x_refresh (int n, char * buffer, int fx, int fy, int fz)
{

  DEBUG ("CommBlock::x_refresh()");
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();

  if ( n != 0) {

    // n == 0 is the call from self to ensure x_refresh()
    // always gets called at least once

    bool gx,gy,gz;
    gx = false;
    gy = false;
    gz = false;

    FieldFace field_face(field_block(), field_descr);

    field_face.set_face(fx,fy,fz);
    field_face.set_ghost(gx,gy,gz);

    field_face.store (n, buffer);
  }

  //--------------------------------------------------
  // Count incoming faces
  // (SHOULD NOT RECOMPUTE EVERY CALL)
  //--------------------------------------------------

  int nx,ny,nz;
  field_block()->size (&nx,&ny,&nz);

  // Determine axes that may be neighbors

  bool fxm = nx > 1;
  bool fxp = nx > 1;
  bool fym = ny > 1;
  bool fyp = ny > 1;
  bool fzm = nz > 1;
  bool fzp = nz > 1;

  // Adjust for boundary faces

  bool periodic = simulation->problem()->boundary()->is_periodic();

  Hierarchy * hierarchy = simulation->hierarchy();

  double lower[3], upper[3];
  hierarchy->lower(&lower[0],&lower[1],&lower[2]);
  hierarchy->upper(&upper[0],&upper[1],&upper[2]);

  bool is_boundary[3][2];
  is_on_boundary (lower,upper,is_boundary);

  fxm = fxm && (periodic || ! is_boundary[axis_x][face_lower]);
  fxp = fxp && (periodic || ! is_boundary[axis_x][face_upper]);
  fym = fym && (periodic || ! is_boundary[axis_y][face_lower]);
  fyp = fyp && (periodic || ! is_boundary[axis_y][face_upper]);
  fzm = fzm && (periodic || ! is_boundary[axis_z][face_lower]);
  fzp = fzp && (periodic || ! is_boundary[axis_z][face_upper]);

  // Count total expected number of incoming faces

  // self

  int count = 1;

  // faces

  if (field_descr->refresh_face(2)) {
    if ( fxm ) ++count;
    if ( fxp ) ++count;
    if ( fym ) ++count;
    if ( fyp ) ++count;
    if ( fzm ) ++count;
    if ( fzp ) ++count;
  }

  // edges

  if (field_descr->refresh_face(1)) {
    if ( fxm && fym ) ++count;
    if ( fxm && fyp ) ++count;
    if ( fxp && fym ) ++count;
    if ( fxp && fyp ) ++count;
    if ( fym && fzm ) ++count;
    if ( fym && fzp ) ++count;
    if ( fyp && fzm ) ++count;
    if ( fyp && fzp ) ++count;
    if ( fzm && fxm ) ++count;
    if ( fzm && fxp ) ++count;
    if ( fzp && fxm ) ++count;
    if ( fzp && fxp ) ++count;
  }

  // corners

  if (field_descr->refresh_face(0)) {
    if ( fxm && fym && fzm ) ++count;
    if ( fxm && fym && fzp ) ++count;
    if ( fxm && fyp && fzm ) ++count;
    if ( fxm && fyp && fzp ) ++count;
    if ( fxp && fym && fzm ) ++count;
    if ( fxp && fym && fzp ) ++count;
    if ( fxp && fyp && fzm ) ++count;
    if ( fxp && fyp && fzp ) ++count;
  }

  //--------------------------------------------------
  // Compute
  //--------------------------------------------------

  if (++count_refresh_face_ >= count) {
    DEBUG ("CommBlock::x_refresh() calling prepare()");
    count_refresh_face_ = 0;
    prepare();
  } else  DEBUG ("CommBlock::x_refresh() skipping prepare()");
}

//----------------------------------------------------------------------

void CommBlock::p_write (int index_output)
{
  TRACE("OUTPUT CommBlock::p_write()");
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  Output * output = simulation->problem()->output(index_output);

  output->write_block(this,field_descr,0,0,0);

  proxy_patch_.s_write();
}

//======================================================================

#endif /* CONFIG_USE_CHARM */

