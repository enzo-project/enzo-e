// See LICENSE_CELLO file for license and copyright information

/// @file     comm_CommBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-11-18
/// @brief    Communication class associated with Blocks

#include "comm.hpp"

//----------------------------------------------------------------------

CommBlock::CommBlock() throw ()
{
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
void CommBlock::p_compute (int cycle, double time, double dt)
{
  // set_cycle(cycle);
  // set_time(time);
  // set_dt(dt);

  DEBUG3 ("CommBlock::p_compute() cycle %d time %f dt %f",cycle,time,dt);
  compute();
}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

void CommBlock::p_initial()
{
  TRACE("CommBlock::p_initial()");
  Simulation * simulation  = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  // Initialize the block

  allocate(field_descr);

  // Set the Block cycle and time to match Simulation's

  TRACE("CommBlock::p_initial Setting time");
  set_cycle(simulation->cycle());
  set_time (simulation->time());
  set_dt   (simulation->dt());

  // Perform any additional initialization for derived class 

  initialize ();

  // Apply the initial conditions 

  Initial * initial = simulation->problem()->initial();

  initial->enforce_block(this,field_descr, simulation->hierarchy());

  // Continue with Patch::s_initial

  proxy_patch_.s_initial();

}

//----------------------------------------------------------------------

void CommBlock::p_output(CkReductionMsg * msg)
{

  DEBUG("CommBlock::p_output()");
  double * min_reduce = (double * )msg->getData();

  double dt_patch   = min_reduce[0];
  bool   stop_patch = min_reduce[1] == 1.0 ? true : false;

  delete msg;

  set_dt   (dt_patch);

  // WARNING: assumes one patch

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  simulation->update_state(cycle_,time_,dt_patch,stop_patch);
 
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // ??? HOW IS cycle_ and time_ update on all processors ensured before index() calls
  // Simulation::p_output()?  Want last block?
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  

  // "root" block calls Simulation::p_output()
  //  if (index() == 0) {
  //    proxy_simulation.p_output();
  //  }

  // Wait for all blocks to check in before calling Simulation::p_output()
  // for next output

  proxy_patch_.s_output();

}
#endif /* CONFIG_USE_CHARM */

void CommBlock::p_read (int index_initial)
{
  INCOMPLETE("CommBlock::p_read");
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

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
#endif /* CONFIG_USE_CHARM */

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

