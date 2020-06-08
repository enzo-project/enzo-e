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
  cello::refresh(ir_post_)->set_callback(CkIndex_Block::p_compute_continue());
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

  p | schedule_; // pupable
  p | courant_;
  p | ir_post_;
  p | required_fields_; // std::vector<str> required fields

}

//----------------------------------------------------------------------

int Method::add_new_refresh_ ()
{
  // set Method::ir_post_

  const int ghost_depth = 4; // std::max(g3[0],std::max(g3[1],g3[2]));
  const int min_face_rank = 0; // cello::config()->adapt_min_face_rank;

  // Set default refresh object
  Refresh refresh_default
    (ghost_depth,min_face_rank, neighbor_leaf, sync_neighbor, 0);

  return cello::simulation()->new_register_refresh(refresh_default);
}

//----------------------------------------------------------------------

int Method::refresh_id_post() const
{
  return ir_post_;
}

//----------------------------------------------------------------------

void Method::set_schedule (Schedule * schedule) throw()
{
  if (schedule_) delete schedule_;
  schedule_ = schedule;
}

//----------------------------------------------------------------------

void Method::define_fields () throw()
{
  /* Ensure required fields are defined for this method */

  FieldDescr * field_descr = cello::field_descr();

  for (int ifield = 0; ifield < required_fields_.size(); ifield++){
    if( ! field_descr->is_field( required_fields_[ifield] )){
      field_descr->insert_permanent( required_fields_[ifield] );
    }
  }

}

//----------------------------------------------------------------------

void Method::define_group_fields (std::vector<std::string> group_fields,
                                  std::string groupname) throw()
{
  /* Ensure fields are grouped correctly */

  FieldDescr * field_descr = cello::field_descr();

  for (int ifield = 0; ifield < group_fields.size(); ifield++){

    // Maybe just throw error here to keep this fully separate from above
    if( ! field_descr->is_field( required_fields_[ifield] )){
      field_descr->insert_permanent( required_fields_[ifield] );
    }

    if (!(field_descr->groups()->is_in( group_fields[ifield], groupname)) ){
      field_descr->groups()->add( group_fields[ifield], groupname);
    }

  }


}

//======================================================================
