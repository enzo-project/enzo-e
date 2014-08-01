// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-05-02
/// @brief    Implentation of the Method class

#include "simulation.hpp"

//----------------------------------------------------------------------

void Method::initialize_ ( CommBlock * comm_block) throw()
{

  Block * block = comm_block -> block();

  block->lower(&xbm_,&ybm_,&zbm_);
  block->upper(&xbp_,&ybp_,&zbp_);

  Hierarchy * hierarchy = comm_block->simulation()->hierarchy();

  hierarchy->lower(&xdm_,&ydm_,&zdm_);
  hierarchy->upper(&xdp_,&ydp_,&zdp_);

  FieldBlock       * field_block = block      -> field_block();
  const FieldDescr * field_descr = comm_block -> field_descr();

  field_count_ = field_descr->field_count();

  field_name_     .resize(field_count_);
  field_array_    .resize(field_count_);
  field_precision_.resize(field_count_);
  gx_             .resize(field_count_);
  gy_             .resize(field_count_);
  gz_             .resize(field_count_);
  mx_             .resize(field_count_);
  my_             .resize(field_count_);
  mz_             .resize(field_count_);

  field_block -> size(&nbx_,&nby_,&nbz_);

  hierarchy->root_size(&ndx_,&ndy_,&ndz_);

  for (int id=0; id<field_count_; id++) {
    field_descr->ghosts (id,&gx_[id],&gy_[id],&gz_[id]);
    field_name_[id] = field_descr->field_name(id);
    field_id_[field_name_[id]] = id;
    field_array_[id] = field_block->values(id);
    field_precision_[id] = field_descr->precision (id);
    mx_[id] = nbx_+2*gx_[id];
    my_[id] = nby_+2*gy_[id];
    mz_[id] = nbz_+2*gz_[id];
  }


  rank_ = comm_block->rank();

  dt_ = comm_block->dt();

  comm_block->cell_width (&hx_,&hy_,&hz_);

}
//======================================================================

