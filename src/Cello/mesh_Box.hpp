// See LICENSE_CELLO file for license and copyright information

/// @file     Mesh_Box.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-09-16
/// @brief    [\ref Mesh] Declaration of the Box class

#ifndef MESH_BOX_HPP
#define MESH_BOX_HPP

enum BoxType { BoxType_receive = 0, BoxType_extra = 1, BoxType_ignored = -1 };
enum class BlockType { undefined, send, receive, extra, receive_coarse, extra_coarse, none };

class Box {

  /// @class    Box
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh]

public: // interface

  /// Constructor
  Box(int rank, int block_size[3], int ghost_depth[3])
    : rank_(rank),
      block_size_(),
      ghost_depth_(),
      ghost_depth_recv_(),
      ghost_depth_send_(),
      coarse_size_(),
      coarse_ghost_(),
      pad_(),
      face_(),
      child_(),
      block_start_(),
      region_start_(),
      region_stop_(),
      centering_()
  {
    for (int i=0; i<2; i++) {
      level_[i] = 0;
    }
    for (int i=0; i<rank_; i++) {
      block_size_[i] = block_size[i];
      ghost_depth_[i] = ghost_depth[i];
      ghost_depth_recv_[i] = ghost_depth[i];
      ghost_depth_send_[i] = 0;
      coarse_size_[i] = block_size[i] / 2;
      coarse_ghost_[i] = ghost_depth_[i]/2 + ghost_depth_[i]%1 + 1;
    }
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    p | rank_;
    PUParray(p,block_size_,3);
    PUParray(p,ghost_depth_,3);
    PUParray(p,ghost_depth_recv_,3);
    PUParray(p,ghost_depth_send_,3);
    PUParray(p,coarse_size_,3);
    PUParray(p,coarse_ghost_,3);
    PUParray(p,level_,2);
    for (int i=0; i<2; i++) {
      PUParray(p,face_[i],3);
      PUParray(p,child_[i],3);
      PUParray(p,block_start_[i],3);
    }
    p | pad_;
    PUParray(p,region_start_,3);
    PUParray(p,region_stop_,3);
    PUParray(p,centering_,3);

  }

  /// Set Block face
  inline void set_block (BoxType bt, int level, int face[3], int child[3])
  {
    level_[bt] = level;
    face_[bt][0] = face[0];
    face_[bt][1] = face[1];
    face_[bt][2] = face[2];
    child_[bt][0] = child[0];
    child_[bt][1] = child[1];
    child_[bt][2] = child[2];
    compute_block_start(bt);
  }

  /// Set number of ghost zones 0 < g <= ghost_depth to receive
  /// Default ghost_depth
  inline void set_recv_ghosts (int ghost_depth_recv[3])
  {
    ghost_depth_recv_[0] = ghost_depth_recv[0];
    ghost_depth_recv_[1] = ghost_depth_recv[1];
    ghost_depth_recv_[2] = ghost_depth_recv[2];
  }

  /// Set number of ghost zones 0 <= g <= ghost_depth to include from the
  /// send block, e.g. for accumulate.  Default 0.
  inline void set_send_ghosts (int ghost_depth_send[3])
  {
    ghost_depth_send_[0] = ghost_depth_send[0];
    ghost_depth_send_[1] = ghost_depth_send[1];
    ghost_depth_send_[2] = ghost_depth_send[2];
  }

  /// Set coarse-zone padding around intersected region when
  /// level = +1
  inline void set_padding (int pad)
  {
    pad_ = pad;
  }

  /// Set coarse field size
  inline void set_coarse (int coarse_size[3], int coarse_ghost[3])
  {
    for (int i=0; i<rank_; i++) {
      coarse_size_[i] = coarse_size[i];
      coarse_ghost_[i] = coarse_ghost[i];
    }
  }

  /// Set field centering
  inline void set_centering (int centering[3])
  {
    for (int i=0; i<rank_; i++) {
      centering_[i] = centering[i];
    }
  }

  /// Determine the start of the specified block
  void compute_block_start(BoxType bt);

  /// Determine send/recv intersection region
  void compute_region();

  /// Get send/box (recv/extra) intersection loop limits for the
  /// specified block
  bool get_start_stop
  (int index_start[3], int index_stop[3],
   BlockType block_intersect, BlockType block_coords, bool lpad);

  bool get_start_size
  (int index_start[3], int index_size[3],
   BlockType block_intersect, BlockType block_coords, bool lpad);

  void get_region_size (int ma3[3])
  {
    for (int i=0; i<rank_; i++) {
      ma3[i] = region_stop_[i] - region_start_[i];
    }
    for (int i=rank_; i<3; i++) {
      ma3[i] = 1;
    }
  }

  void print(const char * mesg)
  {
    CkPrintf ("BOX------ %s -------\n",mesg);
    CkPrintf ("BOX: rank %d\n",rank_);
    CkPrintf ("BOX: block_size_ %d %d %d\n",
              block_size_[0],block_size_[1],block_size_[2]);
    CkPrintf ("BOX: ghost_depth_ %d %d %d\n",
              ghost_depth_[0],ghost_depth_[1],ghost_depth_[2]);
    CkPrintf ("BOX: ghost_depth_recv_ %d %d %d\n",
              ghost_depth_recv_[0],ghost_depth_recv_[1],ghost_depth_recv_[2]);
    CkPrintf ("BOX: ghost_depth_send_ %d %d %d\n",
              ghost_depth_send_[0],ghost_depth_send_[1],ghost_depth_send_[2]);
    CkPrintf ("BOX: pad_ %d\n", pad_);
    CkPrintf ("BOX\n");
    for (int i=0; i<2; i++) {
      CkPrintf ("BOX: level_[%d] %d\n", i,level_[i]);
      CkPrintf ("BOX: face_[%d] %d %d %d\n", i,face_[i][0],face_[i][1],face_[i][2]);
      CkPrintf ("BOX: child_[%d] %d %d %d\n", i,child_[i][0],child_[i][1],child_[i][2]);
      CkPrintf ("BOX: block_start_[%d] %d %d %d\n",
                i,block_start_[i][0],block_start_[i][1],block_start_[i][2]);
    }
    CkPrintf ("BOX\n");
    CkPrintf ("BOX: region_start_ %d %d %d\n",
              region_start_[0],region_start_[1],region_start_[2]);

    CkPrintf ("BOX: region_stop_ %d %d %d\n",
              region_stop_[0],region_stop_[1],region_stop_[2]);
    CkPrintf ("Box----------------\n");
  }

private: // methods

  /// Restrict limits to block, including or excluding ghosts
  bool intersect_regions_
  (int index_min[3], int index_max[3],
   const int index_start[3], const int index_stop[3]);

  /// Apply padding, if any, to the index region
  void apply_padding_(int index_min[3], int index_max[3]);

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Spacial dimension of the boxes, at most three
  int rank_;

  /// Size of blocks
  int block_size_[3];

  /// Depth of allocated ghost zones
  int ghost_depth_[3];

  /// Number of ghost zones receive block wants
  int ghost_depth_recv_[3];

  /// Number of extra send-block ghost zones to include, typically
  /// only used with "accumulate" refresh (default is 0) since otherwise
  /// will overwrite non-ghost zones in receiver
  int ghost_depth_send_[3];

  /// Size of coarse fields
  int coarse_size_[3];

  /// Depth of coarse field ghosts
  int coarse_ghost_[3];

  /// Number of extra coarse zones of padding for interpolation
  /// around the intersected region
  int pad_;

  /// Relative refinement level of recv block relative to the send block
  int level_[2];

  /// Face along which the receive-block is located (-1,0,1)^3
  int face_[2][3];

  /// Child index of send-block or receive block, whichever is finer
  /// (not accessed if both in same level)
  int child_[2][3];

  /// Starting index of the neighbor (send) block relative to the
  /// receive block
  int block_start_[2][3];

  /// Starting and stopping indices of the send-recv intersection region
  int region_start_[3];
  int region_stop_[3];

  /// Centering of field variables in cell (0 = centered; 1 = non-centered)
  int centering_[3];
};

#endif /* MESH_BOX_HPP */

