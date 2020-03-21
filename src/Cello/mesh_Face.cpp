// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Face.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-22
/// @brief    Implementation of the Face class

#include "mesh.hpp"

// #define DEBUG_FACE

//----------------------------------------------------------------------

void Face::adjacency (bool *lx, bool *ly, bool *lz) const
{
  int a1x,a1y,a1z;
  int a2x,a2y,a2z;
  int t1x,t1y,t1z;
  int t2x,t2y,t2z;

  // get mesh levels of block and neighbor indices

  if (index_block_ != index_neighbor_) {
    const int level_1 = index_block_.level();
    const int level_2 = index_neighbor_.level();

    ASSERT2 ("Face::adjacency()",
             "level_1 %d must be within one level of level_2 %d",
             level_1,level_2,
             (std::abs(level_1-level_2)<=1));
  
    // extract array indices of block and neighbor Indices

    const int max_level = std::max(level_1,level_2);

    index_block_.array (&a1x,&a1y,&a1z);
    index_block_.tree  (&t1x,&t1y,&t1z,max_level);
    index_neighbor_.array (&a2x,&a2y,&a2z);
    index_neighbor_.tree  (&t2x,&t2y,&t2z,max_level);

#ifdef DEBUG_FACE    
    CkPrintf ("DEBUG block    level %d\n",level_1);
    CkPrintf ("DEBUG block    array %d %d %d\n",a1x,a1y,a1z);
    CkPrintf ("DEBUG block    tree 0x%08x  0x%08x  0x%08x\n",t1x,t1y,t1z);
    CkPrintf ("DEBUG neighbor level %d\n",level_2);
    CkPrintf ("DEBUG neighbor array %d %d %d\n",a2x,a2y,a2z);
    CkPrintf ("DEBUG block    tree 0x%08x  0x%08x  0x%08x\n",t2x,t2y,t2z);
#endif
    
    // left face array+tree bit indices in finest level

    const int b1lx = (a1x*(1<<max_level) + t1x);
    const int b1ly = (a1y*(1<<max_level) + t1y);
    const int b1lz = (a1z*(1<<max_level) + t1z);

    const int b2lx = (a2x*(1<<max_level) + t2x);
    const int b2ly = (a2y*(1<<max_level) + t2y);
    const int b2lz = (a2z*(1<<max_level) + t2z);

    const int r1 = (level_1 == level_2 - 1) ? 2:1;
    const int r2 = (level_2 == level_1 - 1) ? 2:1;

    // right face array+tree bit indices in finest level

    const int b1rx = b1lx + r1;
    const int b1ry = b1ly + r1;
    const int b1rz = b1lz + r1;

    const int b2rx = b2lx + r2;
    const int b2ry = b2ly + r2;
    const int b2rz = b2lz + r2;

    // Whether block extents are disjoint along each axis
    bool djx = (b1rx < b2lx) || (b2rx < b1lx);
    bool djy = (b1ry < b2ly) || (b2ry < b1ly);
    bool djz = (b1rz < b2lz) || (b2rz < b1lz);

    // Adjust whether disjoint if adjacent on periodic domain
    if (px_) {
      const int axs = ax_ << max_level;
      djx = djx && ((axs + b1lx - b2rx) % axs != 0)
        &&         ((axs + b2lx - b1rx) % axs != 0);
    }
    if (py_) {
      const int ays = ay_ << max_level;
      djy = djy && ((ays + b1ly - b2ry) % ays != 0)
        &&         ((ays + b2ly - b1ry) % ays != 0);
    }
    if (pz_) {
      const int azs = az_ << max_level;
      djz = djz && ((azs + b1lz - b2rz) % azs != 0)
        &&         ((azs + b2lz - b1rz) % azs != 0);
    }

    // Whether a block's extents is contained within another along each axis
    const bool inx = ((b1lx <= b2lx) && (b2rx <= b1rx))
      ||             ((b2lx <= b1lx) && (b1rx <= b2rx));
    const bool iny = ((b1ly <= b2ly) && (b2ry <= b1ry))
      ||             ((b2ly <= b1ly) && (b1ry <= b2ry));
    const bool inz = ((b1lz <= b2lz) && (b2rz <= b1rz))
      ||             ((b2lz <= b1lz) && (b1rz <= b2rz));

#ifdef DEBUG_FACE  
    CkPrintf ("DEBUG disjoint = %d %d %d\n",djx?1:0,djy?1:0,djz?1:0);
    CkPrintf ("DEBUG interior = %d %d %d\n",inx?1:0,iny?1:0,inz?1:0);
#endif  

    if (lx) (*lx) = inx;
    if (ly) (*ly) = iny;
    if (lz) (*lz) = inz;

  } else { // index_block == index_neighbor

    ASSERT1 ("Face::adjacency()", "axis = %d", axis_,
             (0 <= axis_ && axis_ < rank_));
    ASSERT1 ("Face::adjacency()", "face = %d", face_,
             (face_ == -1 || face_ == +1));
    if (lx) (*lx) = (axis_ != 0);
    if (ly) (*ly) = (axis_ != 1);
    if (lz) (*lz) = (axis_ != 2);
    
  }
}

//----------------------------------------------------------------------

void Face::get_offset (int *p_ix, int *p_iy, int *p_iz,
                      int bx, int by, int bz)
{
  // Set default starting index (0,0,0)
  if (p_ix) (*p_ix) = 0;
  if (p_iy) (*p_iy) = 0;
  if (p_iz) (*p_iz) = 0;

  const int level_1 = index_block_.level();
  const int level_2 = index_neighbor_.level();

#ifdef DEBUG_FACE    
  CkPrintf ("DEBUG_FACE_OFFSET levels %d %d\n",level_1,level_2);
#endif
  if (level_1 < level_2) {
    const int max_level = std::max(level_1,level_2);
    int tx,ty,tz;
    index_neighbor_.tree  (&tx,&ty,&tz,max_level);
    // If block is coarse with fine neighbor, return offset of fine
    // face in coarse face using tree bits
#ifdef DEBUG_FACE    
    CkPrintf ("DEBUG_FACE_OFFSET axis %d tree %d %d %d\n",axis_,tx,ty,tz);
#endif
    if (p_ix) (*p_ix) = (axis_ != 0) ? (tx & 1)*bx/2 : 0;
    if (p_iy) (*p_iy) = (axis_ != 1 && rank_ >= 2) ? (ty & 1)*by/2 : 0;
    if (p_iz) (*p_iz) = (axis_ != 2 && rank_ >= 3) ? (tz & 1)*bz/2 : 0;
  }

}
