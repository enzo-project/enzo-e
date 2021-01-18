// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-09-21
/// @brief    

#include "mesh.hpp"

// #define DEBUG_BOX

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
#ifdef DEBUG_BOX  
  CkPrintf ("DEBUG_BOX block start block_size %d %d %d\n",
            block_size_[0],block_size_[1],block_size_[2]);
  CkPrintf ("DEBUG_BOX child[%d]              %d %d %d\n",
            bt,child_[bt][0],child_[bt][1],child_[bt][2]);
  CkPrintf ("DEBUG_BOX face[%d]               %d %d %d\n",
            bt,face_[bt][0],face_[bt][1],face_[bt][2]);
  CkPrintf ("DEBUG_BOX block_start[%d]        %d %d %d\n",
            bt,block_start_[bt][0],block_start_[bt][1],block_start_[bt][2]);
#endif  
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
#ifdef DEBUG_BOX
  CkPrintf ("DEBUG_BOX compute region block_size %d %d %d\n",
            block_size_[0],block_size_[1],block_size_[2]);
  CkPrintf ("DEBUG_BOX compute region child %d %d %d\n",
            child_[0][0],child_[0][1],child_[0][2]);
  CkPrintf ("DEBUG_BOX compute region face %d %d %d\n",
            face_[0][0],face_[0][1],face_[0][2]);
  CkPrintf ("DEBUG_BOX compute region L %d pad %d\n", level_[0],pad_);  
            

  CkPrintf ("DEBUG_BOX compute region region_start-pre %d %d %d\n",
            region_start_[0],region_start_[1],region_start_[2]);
  CkPrintf ("DEBUG_BOX compute region region_stop %d %d %d\n",
            region_stop_[0],region_stop_[1],region_stop_[2]);
#endif  
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
#ifdef DEBUG_BOX
  CkPrintf ("DEBUG_BOX compute region region_start-post %d %d %d\n",
            region_start_[0],region_start_[1],region_start_[2]);
  CkPrintf ("DEBUG_BOX compute region region_stop %d %d %d\n",
            region_stop_[0],region_stop_[1],region_stop_[2]);
#endif  
}

//----------------------------------------------------------------------

bool Box::get_limits
( BoxType bt,  BlockType block_type,  int region_start[3], int region_stop[3])
{
  if (block_type == BlockType::send) {

    for (int i=0; i<rank_; i++) {
      region_start[i] = ghost_depth_[i] + region_start_[i];
      region_stop[i]  = ghost_depth_[i] + region_stop_[i];
    }

  } else if (block_type == BlockType::receive ||
             block_type == BlockType::extra) {

    const float r = std::pow(2.0,1.0*level_[bt]);
  
    for (int i=0; i<rank_; i++) {
      region_start[i] = ghost_depth_[i] +
        r*(region_start_[i] - block_start_[bt][i]);
      region_stop[i]  = ghost_depth_[i] +
        r*(region_stop_[i] - block_start_[bt][i]);
    }

  } else if (block_type == BlockType::array) {

    // receiving coarse-level array for interpolation

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
#ifdef DEBUG_BOX
  CkPrintf ("DEBUG_BOX get_limits type %d  im3 %d %d %d  ip3 %d %d %d\n",
            block_type,
            region_start[0],region_start[1],region_start[2],
            region_stop[0],region_stop[1],region_stop[2]);
#endif  

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

