// See LICENSE_CELLO file for license and copyright information

/// @file     Mesh_Box.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-09-16
/// @brief    [\ref Mesh] Declaration of the Box class

#ifndef MESH_BOX_HPP
#define MESH_BOX_HPP

class Box {

  /// @class    Box
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 
  
public: // interface

  /// Constructor
  Box() throw()
  : rank_(0),
    level_(0),
    r_(0),
    nx_(0), ny_(0), nz_(0),
    gx_(0), gy_(0), gz_(0),
    ix0_(0), iy0_(0), iz0_(0),
    ixm_(0),iym_(0),izm_(0),
    ixp_(0),iyp_(0),izp_(0),
    fx_(0), fy_(0), fz_(0),
    cx_(0), cy_(0), cz_(0),
    gxr_(0), gyr_(0), gzr_(0),
    gxs_(0), gys_(0), gzs_(0)
  {
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    p | rank_;
    p | level_;
    
    p | r_;
    p | nx_;    p | ny_;    p | nz_;
    p | gx_;    p | gy_;    p | gz_;

    p | ix0_;    p | iy0_;    p | iz0_;
    p | ixm_;    p | iym_;    p | izm_;
    p | ixp_;    p | iyp_;    p | izp_;
    p | fx_;    p | fy_;    p | fz_;
    p | cx_;    p | cy_;    p | cz_;
    p | gxr_;   p | gyr_;   p | gzr_;
    p | gxs_;   p | gys_;   p | gzs_;

  }

  /// Set the rank of the problem
  inline void set_rank ( int rank)
  { rank_ = rank; }
  
  /// Set the /relative/ level of the receive block
  inline void set_level (int level)
  {
    level_ = level;

    float r;
    if (level==-1) {
      r = 2.0;
    } else if (level==0) {
      r = 1.0;
    }  else {
      r = 0.5;
    }
    
    r_ = r;
  }
  
  /// Set the face adjacent to the receive-block
  inline void set_face (int fx, int fy, int fz)
  {
    fx_ = fx;
    fy_ = fy;
    fz_ = fz;
  }

  /// Set the child associated with the send-block (if level==+1) or
  /// recv-block (if level==-1)
  inline void set_child (int cx, int cy, int cz)
  {
    cx_ = cx;
    cy_ = cy;
    cz_ = cz;
  }

  /// Set size of blocks
  inline void set_block_size (int nx, int ny, int nz)
  {
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
  }

  /// Set allocated ghost zone depths
  inline void set_ghost_depth (int gx, int gy, int gz)
  {
    gx_ = gx;
    gy_ = gy;
    gz_ = gz;
    set_recv_ghosts(gx,gy,gz);
  }

  /// Set number of ghost zones 0 < g <= ghost_depth to receive
  /// Default ghost_depth
  inline void set_recv_ghosts (int gx, int gy, int gz)
  {
    gxr_ = gx;
    gyr_ = gy;
    gzr_ = gz;
  }
  
  /// Set number of ghost zones 0 <= g <= ghost_depth to include from the
  /// send block, e.g. for accumulate.  Default 0.
  inline void set_send_ghosts (int gx, int gy, int gz)
  {
    gxs_ = gx;
    gys_ = gy;
    gzs_ = gz;
  }

  /// Get recv-send intersection loop limits for send-block
  void get_send_limits (int * ixm, int * ixp,
                        int * iym, int * iyp,
                        int * izm, int * izp);
  
  /// Get recv-send intersection loop limits for recv-block
  void get_recv_limits (int * ixm, int * ixp,
                        int * iym, int * iyp,
                        int * izm, int * izp);

  /// Determine intersection region
  void compute_region();


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Spacial dimension of the boxes, at most three
  int rank_;
  
  /// Relative refinement level of recv block relative to the send block
  int level_;

  /// Refinement ratio between recv- and send-blocks 0.5, 1.0, or 2.0
  float r_;

  /// Size of blocks
  int nx_, ny_, nz_;

  /// Depth of allocated ghost zones
  int gx_, gy_, gz_;

  /// Starting index of the neighbor (send) block relative to the
  /// receive block
  int ix0_, iy0_, iz0_;
  
  /// Starting and stopping indices of the intersected region
  int ixm_, iym_, izm_;
  int ixp_, iyp_, izp_;
  
  /// Face along which the receive-block is located (-1,0,1)^3
  int fx_, fy_, fz_;

  /// Child index of send-block or receive block, whichever is finer
  /// (not accessed if both in same level)
  int cx_, cy_, cz_;

  /// Number of ghost zones receive block wants
  int gxr_, gyr_, gzr_;

  /// Number of extra send-block ghost zones to include, typically
  /// only used with "accumulate" refresh (default is 0) since otherwise
  /// will overwrite non-ghost zones in receiver
  int gxs_, gys_, gzs_;

};

#endif /* MESH_BOX_HPP */

