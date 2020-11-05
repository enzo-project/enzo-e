// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-09-21
/// @brief    

#include "mesh.hpp"
//----------------------------------------------------------------------

bool Box::get_limits (int * im3, int * ip3, BlockType block_type)
{
  bool l_intersect = true;

  if (block_type == BlockType::send) {

    for (int i=0; i<rank_; i++) {
      im3[i] = im3_[i] + g3_[i];
      ip3[i] = ip3_[i] + g3_[i];
    }
    for (int i=rank_; i<3; i++) {
      im3[i] = 0;
      ip3[i] = 1;
    }

  } else if (block_type == BlockType::receive) {

    const float r = std::pow(2.0,1.0*level_);
  
    for (int i=0; i<rank_; i++) {
      im3[i] = r*(im3_[i] - i03_[i]) + g3_[i];
      ip3[i] = r*(ip3_[i] - i03_[i]) + g3_[i];
    }
    for (int i=rank_; i<3; i++) {
      im3[i] = 0;
      ip3[i] = 1;
    }

  } else if (block_type == BlockType::extra) {

    const float ri = std::pow(2.0,-1.0*level_);

    for (int i=0; i<rank_; i++) {
      im3[i] = std::max(im3_[i],i03_[i]);
      ip3[i] = std::min(ip3_[i],int(i03_[i]+ri*n3_[i]));
      l_intersect = l_intersect && (im3[i] < ip3[i]); 
    }

  }

  return l_intersect;
}

//----------------------------------------------------------------------
  
void Box::compute_region()
{
  compute_block_start();

  const float ri = std::pow(2.0,-1.0*level_);
  
  for (int i=0; i<rank_; i++) {
    im3_[i] = std::max
      (-gs3_[i], int(floor(i03_[i] - ri*gr3_[i])));
    ip3_[i] = std::min
      (n3_[i] + gs3_[i], int(ceil(i03_[i] + ri*(n3_[i]+gr3_[i]))));
    if (level_ == +1 && pad_ > 0) {
      im3_[i] -= pad_;
      ip3_[i] += pad_;
    }
  }
}

//======================================================================

void Box::compute_block_start()
{
  // compute start i03_[] of the receive (or extra) block given
  // relative level_, f3_[], and c3_[], defining the block, relative
  // to the sending block
  if (level_ == -1) {
    for (int i=0; i<rank_; i++) {
      i03_[i] = n3_[i]*(2.0*f3_[i] - c3_[i]);
    }
  } else if (level_ == 0) {
    for (int i=0; i<rank_; i++) {
      i03_[i] = n3_[i]*f3_[i];
    }
  } else if (level_ == +1) {
    for (int i=0; i<rank_; i++) {
      i03_[i] = n3_[i]*(f3_[i] + 0.5*c3_[i]);
    }
  }
}

