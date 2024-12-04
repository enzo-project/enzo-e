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
    courant_(courant),
    neighbor_type_(neighbor_leaf),
    max_supercycle_(1),
    index_method_(-1),
    super_field_(),
    super_field_curr_(),
    super_field_prev_(),
    is_time_curr_(-1),
    is_time_prev_(-1)
{
  ir_post_ = add_refresh_();
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
  p | neighbor_type_;
  p | max_supercycle_;
  p | index_method_;
  p | super_field_;
  p | super_field_curr_;
  p | super_field_prev_;
  p | is_time_curr_;
  p | is_time_prev_;
}

//----------------------------------------------------------------------

int Method::add_refresh_ (int neighbor_type)
{
  // set Method::ir_post_

  const int ghost_depth = 4; // std::max(g3[0],std::max(g3[1],g3[2]));
  const int min_face_rank = 0; // cello::config()->adapt_min_face_rank;

  // Set default refresh object
  Refresh refresh_default
    (ghost_depth,min_face_rank, neighbor_type, sync_neighbor, 0);

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

//======================================================================

bool Method::is_solve_cycle_(Block * block)
{
  return (block->state()->method(index()).step() == 0);
}

//----------------------------------------------------------------------

int Method::super_define_fields_
(std::string field,
 std::string field_curr,
 std::string field_prev)
{
  // define fields 
  const int id_field      = cello::define_field(field);
  const int id_field_curr = cello::define_field(field_curr);
  const int id_field_prev = cello::define_field(field_prev);

  // save indices
  super_field_.push_back(id_field);
  super_field_curr_.push_back(id_field_curr);
  super_field_prev_.push_back(id_field_prev);

  // save time (if not already saved)
  if (is_time_curr_ == -1) {
    ScalarDescr * scalar_descr_double = cello::scalar_descr_double();
    is_time_curr_ = scalar_descr_double->new_value(name()+"_time_curr");
    is_time_prev_ = scalar_descr_double->new_value(name()+"_time_prev");
  }

  const int id_super = (super_field_.size() - 1);

  return id_super;
}

//----------------------------------------------------------------------

void Method::super_shift_fields_(Block * block )
{
  Field field = block->data()->field();

  // Copy field_curr to field_prev

  for (size_t id_super=0; id_super < super_field_.size(); id_super++) {

    const int id_field_curr = super_field_curr_[id_super];
    const int id_field_prev = super_field_prev_[id_super];
    cello_float * values_curr = (cello_float*) field.values (id_field_curr);
    cello_float * values_prev = (cello_float*) field.values (id_field_prev);
    ASSERT1("EnzoMethodGravity::super_save_field_()",
           "field_[curr|prev] fields not defined for id_super = %d",
            id_super,
           (values_curr && values_prev));

    int mx,my,mz;
    field.dimensions (id_field_curr,&mx,&my,&mz);
    const int m = mx*my*mz;
    for (int i=0; i<m; i++) {
      values_prev[i] = values_curr[i];
    }
  }

  // Copy time_curr to time_prev

  double &     time_prev = *block->data()->scalar_double().value(is_time_prev_);
  const double time_curr = *block->data()->scalar_double().value(is_time_curr_);
  time_prev = time_curr;
}

//----------------------------------------------------------------------

void Method::super_save_fields_(Block * block )
{
  Field field = block->data()->field();

  // Copy field_curr to field_prev

  for (size_t id_super=0; id_super < super_field_.size(); id_super++) {

    const int id_field_curr = super_field_curr_[id_super];
    const int id_field      = super_field_[id_super];
    cello_float * values_curr = (cello_float*) field.values (id_field_curr);
    cello_float * values      = (cello_float*) field.values (id_field);
    ASSERT1("EnzoMethodGravity::super_shift_field_()",
           "field_[curr|prev] fields not defined for id_super = %d",
            id_super,
           (values_curr && values));

    int mx,my,mz;
    field.dimensions (id_field_curr,&mx,&my,&mz);
    const int m = mx*my*mz;
    for (int i=0; i<m; i++) {
      values_curr[i] = values[i];
    }
  }

  // Copy time to time_curr

  double &     time_curr = *block->data()->scalar_double().value(is_time_curr_);
  const double time      = block->state()->time();
  time_curr = time;
}

//----------------------------------------------------------------------

void Method::super_extrapolate_fields_(Block * block, double time )
{
  // Get current and previous times

  auto scalar_double(block->data()->scalar_double());
  const double tc = *scalar_double.value(is_time_curr_);
  const double tp = *scalar_double.value(is_time_prev_);
  const double t  = time;

  // compute extrapolation coefficients

  const double cp = (tc - t) / (tc - tp);
  const double cc = (t - tp) / (tc - tp);

  // perform potential extrapolation on all fields

  for (size_t id_super=0; id_super < super_field_.size(); id_super++) {

    const int id_field      = super_field_[id_super];
    const int id_field_curr = super_field_curr_[id_super];
    const int id_field_prev = super_field_prev_[id_super];

    Field field = block->data()->field();
    cello_float * values      = (cello_float*) field.values (id_field);
    cello_float * values_curr = (cello_float*) field.values (id_field_curr);
    cello_float * values_prev = (cello_float*) field.values (id_field_prev);

    ASSERT1("Method::super_extrapolate_fields_()",
            "Missing field %d", id_field, values);
    ASSERT1("Method::super_extrapolate_fields_()",
            "Missing field_curr %d", id_field_curr, values_curr);
    ASSERT1("Method::super_extrapolate_fields_()",
            "Missing field_prev %d", id_field_prev, values_prev);

    const int m = field.dimensions (id_field);
    for (int i=0; i<m; i++) {
      values[i] = cc*values_curr[i] + cp*values_prev[i];
    }
  }
}

//----------------------------------------------------------------------

void Method::super_update_time_(Block * block, double time)
{
  auto scalar_double(block->data()->scalar_double());
  *scalar_double.value(is_time_curr_) = time;
}
