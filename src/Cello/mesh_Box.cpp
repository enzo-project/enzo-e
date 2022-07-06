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
}

//----------------------------------------------------------------------

void Box::compute_region()
{
  const float ri = std::pow(2.0,-1.0*level_[0]);

  // initialize default for rank < 3
  for (int i=0; i<3; i++) {
    region_start_[i] = 0;
    region_stop_[i] = 1;
  }
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
}

//----------------------------------------------------------------------

bool Box::get_start_size
( int index_start[3], int index_size[3],
  BlockType block_intersect, BlockType block_coords, bool lpad)
{
  bool intersect = get_start_stop
    (index_start,index_size, block_intersect,block_coords,lpad);
  for (int i=0; i<rank_; i++) {
    index_size[i] -= index_start[i];
  }
  return intersect;
}

//----------------------------------------------------------------------

bool Box::get_start_stop
( int index_min[3], int index_max[3],
  BlockType block_intersect, BlockType block_coords,
  bool lpad)
{
  if (pad_ != 0) {
    ASSERT ("Box::get_start_stop",
            "Expecting lpad = false with block_coords == receive",
            ! ((lpad == true) && (block_coords == BlockType::receive)));
    ASSERT ("Box::get_start_stop",
            "Expecting lpad = true with block_coords == send",
            ! ((lpad == false) && (block_coords == BlockType::send)));
    ASSERT ("Box::get_start_stop",
            "Expecting lpad = true with block_coords == *_coarse",
            ! ((lpad == false) &&
             ((block_coords == BlockType::receive_coarse) ||
              (block_coords == BlockType::extra_coarse))));
  }
  // initialize index_[min|max] to be the S->r region
  for (int i=0; i<3; i++) {
    index_min[i] = region_start_[i];
    index_max[i] = region_stop_[i];
  }

  // apply padding if needed
  if (lpad) apply_padding_(index_min,index_max);

  int n3[3] = {1,1,1};
  int g3[3] = {0,0,0};
  for (int i=0; i<rank_; i++) {
    n3[i] = block_size_[i];
    g3[i] = ghost_depth_[i];
  }

  // define intersecting block limits
  int block_min[3] = {0,0,0};
  int block_max[3] = {1,1,1};

  if (block_intersect == BlockType::none) {

    for (int i=0; i<rank_; i++) {
      block_min[i] = -std::numeric_limits<int>::max();
      block_max[i] = +std::numeric_limits<int>::max();
    }

  } else if (block_intersect == BlockType::send) {

    for (int i=0; i<rank_; i++) {
      block_min[i] = -ghost_depth_send_[i];
      block_max[i] = n3[i] + ghost_depth_send_[i];
    }

  } else if (block_intersect == BlockType::receive) {

    const float r = std::pow(2.0,-1.0*level_[0]);

    for (int i=0; i<rank_; i++) {
      block_min[i] = block_start_[0][i] + r*(-ghost_depth_recv_[i]);
      block_max[i] = block_start_[0][i] + r*(n3[i]+ghost_depth_recv_[i]);
    }

  } else if (block_intersect == BlockType::extra) {

    const float r = std::pow(2.0,-1.0*level_[1]);

    for (int i=0; i<rank_; i++) {
      block_min[i] = block_start_[1][i] + r*(-ghost_depth_send_[i]);
      block_max[i] = block_start_[1][i] + r*(n3[i]+ghost_depth_send_[i]);
    }

  } else if (block_intersect == BlockType::receive_coarse ||
             block_intersect == BlockType::extra_coarse) {

    ERROR ("Box::get_start_stop",
           "Intersecting block cannot be the coarse array");

  }

  // compute the block - region intersection
  const bool l_intersect = intersect_regions_
    (index_min,index_max, block_min,block_max);

  // transform to block_coords
  if (block_coords == BlockType::send) {

    for (int i=0; i<rank_; i++) {
      index_min[i] += g3[i];
      index_max[i] += g3[i];
    }

  } else if (block_coords == BlockType::receive) {

    const float r = std::pow(2.0,1.0*level_[0]);
    for (int i=0; i<rank_; i++) {
      index_min[i] = g3[i] + r*(index_min[i] - block_start_[0][i]);
      index_max[i] = g3[i] + r*(index_max[i] - block_start_[0][i]);
    }

  } else if (block_coords == BlockType::extra) {

    const float r = std::pow(2.0,1.0*level_[1]);
    for (int i=0; i<rank_; i++) {
      index_min[i] = g3[i] + r*(index_min[i] - block_start_[1][i]);
      index_max[i] = g3[i] + r*(index_max[i] - block_start_[1][i]);
    }

  } else if (block_coords == BlockType::receive_coarse) {

    for (int i=0; i<rank_; i++) {
      index_min[i] = coarse_ghost_[i] + (index_min[i] - block_start_[0][i]);
      index_max[i] = coarse_ghost_[i] + (index_max[i] - block_start_[0][i]);
    }

  } else if (block_coords == BlockType::extra_coarse) {

    for (int i=0; i<rank_; i++) {
      index_min[i] = coarse_ghost_[i] + (index_min[i] - block_start_[1][i]);
      index_max[i] = coarse_ghost_[i] + (index_max[i] - block_start_[1][i]);
    }

  }

  for (int i=0; i<rank_; i++) {
    if (centering_[i]) {
      if (face_[0][i] == 0) {
        ++index_max[i];
      } else if (face_[0][i] == -1) {
        ++index_min[i];
        ++index_max[i];
      }
    }
  }
  // return whether region intersects
  return l_intersect;
}

//----------------------------------------------------------------------

bool Box::intersect_regions_
(int index_min[3], int index_max[3],
 const int block_min[3], int const block_max[3])
{

  bool l_intersect = true;
  for (int i=0; i<rank_; i++) {
    index_min[i] = std::max(index_min[i],block_min[i]);
    index_max[i] = std::min(index_max[i], block_max[i]);
    l_intersect = l_intersect && (index_min[i] < index_max[i]);
  }
  return l_intersect;
}


void Box::apply_padding_(int index_start[3], int index_stop[3])
{
  for (int i=0; i<rank_; i++) {
    index_start[i] -= pad_;
    index_stop[i]  += pad_;
  }
}
