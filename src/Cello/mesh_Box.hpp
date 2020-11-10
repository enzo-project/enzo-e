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

  enum class BlockType { send, receive, extra, array };
  
  /// Constructor
  Box() throw()
  : rank_(0),
    n3_(),
    g3_(),
    f3_(),
    c3_(),
    level_(0),
    i03_(),
    im3_(),
    ip3_(),
    gr3_(),
    gs3_(),
    pad_(0)
  {
  }

  /// Constructor
  Box(int rank, int n3[3], int g3[3])
    : rank_(rank),
      n3_(),
      g3_(),
      f3_(),
      c3_(),
      level_(0),
      i03_(),
      im3_(),
      ip3_(),
      gr3_(),
      gs3_(),
      pad_()
  {
    for (int i=0; i<rank_; i++) {
      n3_[i] = n3[i];
      g3_[i] = g3[i];
      gr3_[i] = g3[i];
    }
  }
  
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    p | rank_;
    PUParray(p,n3_,3);
    PUParray(p,g3_,3);
    PUParray(p,f3_,3);
    PUParray(p,c3_,3);
    p | level_;
    PUParray(p,i03_,3);
    PUParray(p,im3_,3);
    PUParray(p,ip3_,3);
    PUParray(p,gr3_,3);
    PUParray(p,gs3_,3);
    p | pad_;

  }

  /// Set the rank of the problem
  inline void set_rank ( int rank)
  { rank_ = rank; }
  
  /// Set the /relative/ level of the receive block
  inline void set_level (int level)
  {
    level_ = level;
  }
  
  /// Set the face adjacent to the receive-block
  inline void set_face (int f3[3])
  {
    f3_[0] = f3[0];
    f3_[1] = f3[1];
    f3_[2] = f3[2];
  }

  /// Set the child associated with the send-block (if level==+1) or
  /// recv-block (if level==-1)
  inline void set_child (int c3[3])
  {
    c3_[0] = c3[0];
    c3_[1] = c3[1];
    c3_[2] = c3[2];
  }

  /// Set Block face
  inline void set_block (int level, int f3[3], int c3[3])
  {
    level_ = level;
    f3_[0] = f3[0];
    f3_[1] = f3[1];
    f3_[2] = f3[2];
    c3_[0] = c3[0];
    c3_[1] = c3[1];
    c3_[2] = c3[2];
    compute_block_start();
  }
  
  /// Set size of blocks
  inline void set_block_size (int n3[3])
  {
    n3_[0] = n3[0];
    n3_[1] = n3[1];
    n3_[2] = n3[2];
  }

  /// Set allocated ghost zone depths
  inline void set_ghost_depth (int g3[3])
  {
    g3_[0] = g3[0];
    g3_[1] = g3[1];
    g3_[2] = g3[2];
    set_recv_ghosts(g3);
  }

  /// Set number of ghost zones 0 < g <= ghost_depth to receive
  /// Default ghost_depth
  inline void set_recv_ghosts (int gr3[3])
  {
    gr3_[0] = gr3[0];
    gr3_[1] = gr3[1];
    gr3_[2] = gr3[2];
  }
  
  /// Set number of ghost zones 0 <= g <= ghost_depth to include from the
  /// send block, e.g. for accumulate.  Default 0.
  inline void set_send_ghosts (int gs3[3])
  {
    gs3_[0] = gs3[0];
    gs3_[1] = gs3[1];
    gs3_[2] = gs3[2];
  }

  /// Set coarse-zone padding around intersected region when
  /// level = +1
  inline void set_padding (int pad)
  {
    pad_ = pad;
  }
  
  /// Determine the start of the currently defined receive or extra
  /// block
  void compute_block_start();

  /// Determine intersection region
  void compute_region();

  /// Get recv-send intersection loop limits for send-block
  bool get_limits (int im3[3], int ip3[3], BlockType block_type);

  /// Restrict limits to block, including or excluding ghosts
  void restrict_limits (int im3[3], int ip3[3], bool include_ghosts);

  void print(const char * mesg)
  {
    CkPrintf ("BOX------ %s -------\n",mesg);
    CkPrintf ("BOX: rank %d\n",rank_);
    CkPrintf ("BOX: n3_ %d %d %d\n", n3_[0],n3_[1],n3_[2]);
    CkPrintf ("BOX: g3_ %d %d %d\n", g3_[0],g3_[1],g3_[2]);
    CkPrintf ("BOX: gr3_ %d %d %d\n", gr3_[0],gr3_[1],gr3_[2]);
    CkPrintf ("BOX: gs3_ %d %d %d\n", gs3_[0],gs3_[1],gs3_[2]);
    CkPrintf ("BOX: pad_ %d\n", pad_);
    CkPrintf ("BOX\n");
    CkPrintf ("BOX: level_ %d\n", level_);
    CkPrintf ("BOX: f3_ %d %d %d\n", f3_[0],f3_[1],f3_[2]);
    CkPrintf ("BOX: c3_ %d %d %d\n", c3_[0],c3_[1],c3_[2]);
    CkPrintf ("BOX: i03_ %d %d %d\n", i03_[0],i03_[1],i03_[2]);
    CkPrintf ("BOX\n");
    CkPrintf ("BOX: im3_ %d %d %d\n", im3_[0],im3_[1],im3_[2]);
    CkPrintf ("BOX: ip3_ %d %d %d\n", ip3_[0],ip3_[1],ip3_[2]);
    CkPrintf ("Box----------------\n");
  }
  
private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Spacial dimension of the boxes, at most three
  int rank_;
  
  /// Size of blocks
  int n3_[3];

  /// Depth of allocated ghost zones
  int g3_[3];

  /// Face along which the receive-block is located (-1,0,1)^3
  int f3_[3];

  /// Child index of send-block or receive block, whichever is finer
  /// (not accessed if both in same level)
  int c3_[3];

  /// Relative refinement level of recv block relative to the send block
  int level_;

  /// Starting index of the neighbor (send) block relative to the
  /// receive block
  int i03_[3];
  
  /// Starting and stopping indices of the intersected region
  int im3_[3];
  int ip3_[3];
  
  /// Number of ghost zones receive block wants
  int gr3_[3];

  /// Number of extra send-block ghost zones to include, typically
  /// only used with "accumulate" refresh (default is 0) since otherwise
  /// will overwrite non-ghost zones in receiver
  int gs3_[3];

  /// Number of extra coarse zones of padding for interpolation
  /// around the intersected region
  int pad_;

};

#endif /* MESH_BOX_HPP */

