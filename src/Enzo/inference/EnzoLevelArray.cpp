// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoLevelArray.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-10-28
/// @brief    Implements the EnzoLevelArray class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoLevelArray::EnzoLevelArray
(std::string field_group,
 int level_base, int level_array, int level_infer,
 int nax, int nay, int naz)
  : CBase_EnzoLevelArray(),
    level_base_(level_base),
    level_array_(level_array),
    level_infer_(level_infer),
    nax_(nax),
    nay_(nay),
    naz_(naz),
    nix_(),
    niy_(),
    niz_(),
    field_group_(field_group),
    num_fields_(cello::field_groups()->size(field_group_)),
    field_values_(),
    volume_ratio_(0.0),
    spheres_()
{
  ASSERT2("EnzoLevelArray::EnzoLevelArray()",
          "level_base %d must be between 0 and level_array %d",
          level_base_,level_array_,
          (0 <= level_base_) && (level_base_ <= level_array_));
  ASSERT2("EnzoLevelArray::EnzoLevelArray()",
          "level_infer %d must be at least level_array %d",
          level_infer_,level_array_,
          (level_array_ <= level_infer_));

  const int r = pow(2,level_infer_ - level_array_);
  int nd3[3] = {
    cello::config()->mesh_root_blocks[0],
    cello::config()->mesh_root_blocks[1],
    cello::config()->mesh_root_blocks[2] };
  int md3[3] = {
    cello::config()->mesh_root_size[0],
    cello::config()->mesh_root_size[1],
    cello::config()->mesh_root_size[2] };
  int nb3[3]= { md3[0]/nd3[0], md3[1]/nd3[1], md3[2]/nd3[2]};

  nix_ = r*nb3[0];
  niy_ = r*nb3[1];
  niz_ = r*nb3[2];
  field_values_.resize(num_fields_);
  for (int i=0; i<num_fields_; i++) {
    field_values_[i].resize(nix_*niy_*niz_);
  }
  proxy_enzo_simulation[0].p_infer_array_created();
}

//----------------------------------------------------------------------

void EnzoLevelArray::pup (PUP::er &p)
{

  TRACEPUP;

  CBase_EnzoLevelArray::pup(p);

  p | level_array_;
  p | level_base_;
  p | level_infer_;
  p | nax_;
  p | nay_;
  p | naz_;
  p | nix_;
  p | niy_;
  p | niz_;
  p | field_group_;
  p | num_fields_;
  p | field_values_;
  p | volume_ratio_;
  p | spheres_;
}

//----------------------------------------------------------------------

EnzoLevelArray::~EnzoLevelArray()
{ }

