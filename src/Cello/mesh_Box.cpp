// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-09-21
/// @brief    

#include "mesh.hpp"

// #define DEBUG_BOX

//----------------------------------------------------------------------

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
  for (int i=rank_; i<3; i++) {
    i03_[i] = 0;
  }
#ifdef DEBUG_BOX  
  CkPrintf ("DEBUG_BOX block start n3 %d %d %d c3 %d %d %d f3 %d %d %d:  i03 %d %d %d\n",
            n3_[0],n3_[1],n3_[2],
            c3_[0],c3_[1],c3_[2],
            f3_[0],f3_[1],f3_[2],
            i03_[0],i03_[1],i03_[2]);
#endif  
}

//----------------------------------------------------------------------
  
void Box::compute_region()
{
  const float ri = std::pow(2.0,-1.0*level_);
  
  for (int i=0; i<rank_; i++) {
    im3_[i] = std::max
      (-gs3_[i], int(floor(i03_[i] - ri*gr3_[i])));
    ip3_[i] = std::min
      (n3_[i] + gs3_[i], int(ceil(i03_[i] + ri*(n3_[i]+gr3_[i]))));
  }
  for (int i=rank_; i<3; i++) {
    im3_[i] = 0;
    ip3_[i] = 1;
  }
#ifdef DEBUG_BOX
  CkPrintf ("DEBUG_BOX compute region  n3 %d %d %d c3 %d %d %d f3 %d %d %d L %d pad %d\n",
            n3_[0],n3_[1],n3_[2],c3_[0],c3_[1],c3_[2],f3_[0],f3_[1],f3_[2],level_,pad_);
  CkPrintf ("DEBUG_BOX compute region im3-pre %d %d %d ip3 %d %d %d\n",
            im3_[0],im3_[1],im3_[2],
            ip3_[0],ip3_[1],ip3_[2]);
#endif  
  for (int i=0; i<rank_; i++) {
    if (pad_ > 0) {
      if (level_ == +1) {
        im3_[i] -= pad_;
        ip3_[i] += pad_;
      }
      if (level_ == -1) {
        im3_[i] -= 2*pad_;
        ip3_[i] += 2*pad_;
      }
    }
  }
#ifdef DEBUG_BOX
  CkPrintf ("DEBUG_BOX compute region im3-post %d %d %d ip3 %d %d %d\n",
            im3_[0],im3_[1],im3_[2],
            ip3_[0],ip3_[1],ip3_[2]);
#endif  
}

//----------------------------------------------------------------------

bool Box::get_limits (int im3[3], int ip3[3], BlockType block_type)
{
  if (block_type == BlockType::send) {

    for (int i=0; i<rank_; i++) {
      im3[i] = im3_[i] + g3_[i];
      ip3[i] = ip3_[i] + g3_[i];
    }
    for (int i=rank_; i<3; i++) {
      im3[i] = 0;
      ip3[i] = 1;
    }

  } else if (block_type == BlockType::receive ||
             block_type == BlockType::extra) {

    const float r = std::pow(2.0,1.0*level_);
  
    for (int i=0; i<rank_; i++) {
      im3[i] = r*(im3_[i] - i03_[i]) + g3_[i];
      ip3[i] = r*(ip3_[i] - i03_[i]) + g3_[i];
    }
    for (int i=rank_; i<3; i++) {
      im3[i] = 0;
      ip3[i] = 1;
    }

    if (block_type == BlockType::extra) {
      restrict_limits (im3,ip3,false);
    }
    
  } else if (block_type == BlockType::array) {

    // receiving coarse-level array for interpolation

    const float ri = std::pow(2.0,-1.0*level_);
    
    for (int i=0; i<rank_; i++) {
      im3[i] = std::max(im3_[i],i03_[i]) - im3_[i];
      ip3[i] = std::min(ip3_[i],int(i03_[i]+ri*n3_[i])) - im3_[i];
    }
    for (int i=rank_; i<3; i++) {
      im3[i] = 0;
      ip3[i] = 1;
    }

  }

  bool l_intersect = true;
  for (int i=0; i<rank_; i++) {
    l_intersect = l_intersect && (im3[i] < ip3[i]); 
  }

  return l_intersect;
}

//----------------------------------------------------------------------

void Box::restrict_limits (int im3[3], int ip3[3], bool include_ghosts)
{
  if (include_ghosts) {
    for (int i=0; i<rank_; i++) {
      im3[i] = std::max(im3[i],0);
      ip3[i] = std::min(ip3[i],n3_[i]+2*g3_[i]);
    }
  } else {
    for (int i=0; i<rank_; i++) {
      im3[i] = std::max(im3[i],g3_[i]);
      ip3[i] = std::min(ip3[i],n3_[i]+g3_[i]);
    }
  }
}

