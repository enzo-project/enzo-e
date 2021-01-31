// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-09-21
/// @brief    

#include "mesh.hpp"

//----------------------------------------------------------------------

void Box::compute_block_start(BoxType bt)
{
  // compute start block_start_[] of the receive (or extra) block given
  // relative level_, face_[], and child_[], defining the block, relative
  // to the sending block
  if (level_[bt] == -1) {
    for (int i=0; i<rank_; i++) {
      block_start_[bt][i] = block_size_[i]*(2.0*face_[bt][i] - child_[bt][i]);
    }
  } else if (level_[bt] == 0) {
    for (int i=0; i<rank_; i++) {
      block_start_[bt][i] = block_size_[i]*face_[bt][i];
    }
  } else if (level_[bt] == +1) {
    for (int i=0; i<rank_; i++) {
      block_start_[bt][i] = block_size_[i]*(face_[bt][i] + 0.5*child_[bt][i]);
    }
  }
  for (int i=rank_; i<3; i++) {
    block_start_[bt][i] = 0;
  }
}

//----------------------------------------------------------------------
  
void Box::compute_region()
{
  const float ri = std::pow(2.0,-1.0*level_[0]);
  
  for (int i=0; i<rank_; i++) {
    const int recv_start =
      int(floor(block_start_[0][i] - ri*ghost_depth_recv_[i]));
    const int recv_stop =
      int(ceil(block_start_[0][i] + ri*(block_size_[i]+ghost_depth_recv_[i])));
    region_start_[i] = std::max
      (-ghost_depth_send_[i], recv_start);
    region_stop_[i] = std::min
      (block_size_[i] + ghost_depth_send_[i], recv_stop);
  }
  for (int i=rank_; i<3; i++) {
    region_start_[i] = 0;
    region_stop_[i] = 1;
  }
  if (pad_ > 0) {
    if (level_[0] == +1) {
      for (int i=0; i<rank_; i++) {
        region_start_[i] -= pad_;
        region_stop_[i]  += pad_;
      }
    } else if (level_[0] == -1) {
      for (int i=0; i<rank_; i++) {
        region_start_[i] -= 2*pad_;
        region_stop_[i]  += 2*pad_;
      }
    }
  }
}

//----------------------------------------------------------------------

bool Box::get_limits
( BlockType block_type,  int region_start[3], int region_stop[3], BoxType bt )
{
  if (block_type == BlockType::send) {

    for (int i=0; i<rank_; i++) {
      region_start[i] = ghost_depth_[i] + region_start_[i];
      region_stop[i]  = ghost_depth_[i] + region_stop_[i];
    }

  } else if (block_type == BlockType::receive ||
             block_type == BlockType::extra) {
    ASSERT("Box::get_limits()",
           "Expected BoxType bt to not be BoxType_ignored",
           (bt != BoxType_ignored));
    const float r = std::pow(2.0,1.0*level_[bt]);
  
    for (int i=0; i<rank_; i++) {
      region_start[i] = ghost_depth_[i] +
        r*(region_start_[i] - block_start_[bt][i]);
      region_stop[i]  = ghost_depth_[i] +
        r*(region_stop_[i] - block_start_[bt][i]);
    }

  } else if (block_type == BlockType::array) {

    // receiving coarse-level array for interpolation

    ASSERT("Box::get_limits()",
           "Expected BoxType bt to not be BoxType_ignored",
           (bt != BoxType_ignored));

    const float ri = std::pow(2.0,-1.0*level_[bt]);
    
    for (int i=0; i<rank_; i++) {
      region_start[i] = -region_start_[i] +
        std::max(region_start_[i],block_start_[bt][i]) ;
      region_stop[i]  = -region_start_[i] +
        std::min(region_stop_[i],int(block_start_[bt][i]+ri*block_size_[i])) ;
    }

  }
  for (int i=rank_; i<3; i++) {
    region_start[i] = 0;
    region_stop[i] = 1;
  }

  if (block_type == BlockType::extra) {
    restrict_limits (region_start,region_stop,false);
  }
    
  bool l_intersect = true;
  for (int i=0; i<rank_; i++) {
    l_intersect = l_intersect && (region_start[i] < region_stop[i]); 
  }

  return l_intersect;
}

//----------------------------------------------------------------------

void Box::restrict_limits
(int region_start[3], int region_stop[3], bool include_ghosts)
{
  if (include_ghosts) {
    for (int i=0; i<rank_; i++) {
      const int block_start = 0;
      const int block_stop  = block_size_[i]+2*ghost_depth_[i];
      region_start[i] = std::max(region_start[i],block_start);
      region_stop[i]  = std::min(region_stop[i], block_stop);
    }
  } else {
    for (int i=0; i<rank_; i++) {
      const int block_start = ghost_depth_[i];
      const int block_stop  = ghost_depth_[i] + block_size_[i];;
      region_start[i] = std::max(region_start[i],block_start);
      region_stop[i]  = std::min(region_stop[i], block_stop);
    }
  }
}

