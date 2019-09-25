// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-09-04
/// @brief    

#include "problem.hpp"

double Method::courant_global = 1.0;

//----------------------------------------------------------------------

Method::Method (double courant) throw()
  : schedule_(NULL),
    courant_(courant)
{
  ir_post_ = add_new_refresh_();
  new_refresh(ir_post_).set_callback(CkIndex_Block::r_compute_continue());
}

//----------------------------------------------------------------------

Method::~Method() throw()
{
  delete schedule_;
}

//----------------------------------------------------------------------

void Method::pup (PUP::er &p)
{ TRACEPUP;
  PUP::able::pup(p);

  bool pk = p.isPacking();
  bool up = p.isUnpacking();

  p | schedule_; // pupable
  p | courant_;
  p | ir_post_;

}

//----------------------------------------------------------------------

int Method::add_new_refresh_ ()
{
  // set Method::ir_post_

  const int * g3 = cello::config()->field_ghost_depth;
  const int ghost_depth = 4; // std::max(g3[0],std::max(g3[1],g3[2]));
  const int min_face_rank = 0; // cello::config()->adapt_min_face_rank;

  // Set default refresh object
  Refresh refresh_default
    (ghost_depth,min_face_rank, neighbor_leaf, sync_neighbor, 0);

  return cello::simulation()->new_register_refresh(refresh_default);
}

//----------------------------------------------------------------------

Refresh & Method::new_refresh(int ir)
{
  return cello::simulation()->new_refresh_list(ir);
}

//----------------------------------------------------------------------

int Method::refresh_post_id() const
{
  return ir_post_;
}

//----------------------------------------------------------------------

Refresh & Method::refresh_post()
{
  return cello::simulation()->new_refresh_list(ir_post_);
}

//----------------------------------------------------------------------

void Method::set_schedule (Schedule * schedule) throw()
{ 
  if (schedule_) delete schedule_;
  schedule_ = schedule;
}

//======================================================================

