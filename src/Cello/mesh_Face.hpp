// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Face.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-08
/// @brief    [\ref Mesh] Declaration of the Face class

#ifndef MESH_FACE_HPP
#define MESH_FACE_HPP

class Face {

  /// @class    Face
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh]

public: // interface

  /// Constructor
  Face(Index index_block, Index index_neighbor,
       int rank,
       int ax, int ay=0, int az=0,
       bool px=false, bool py=false, bool pz=false,
       int fx=0, int fy=0, int fz=0) throw()
    :  index_block_(index_block),
       index_neighbor_(index_neighbor),
       ax_(ax),ay_(ay),az_(az),
       px_(px),py_(py),pz_(pz),
       fx_(fx),fy_(fy),fz_(fz),
       axis_(-1), face_(0),
       rank_(rank)
  {
    ASSERT1 ("Face()",
            "rank out of bounds: 1 <= %d <= 3",
            rank_, (1 <= rank_ && rank_ <= 3));
  }

  /// operator < () is used for storing Face objects in a std::map
  bool operator < (const Face  & face) const
  {
    const Index &i1 = this->index_block_;
    const Index &i2 = this->index_neighbor_;
    const Index &j1 = face.index_block_;
    const Index &j2 = face.index_neighbor_;

    bool lt = false;
    if (i1 < j1) {
      // Block index is less than face's
      lt = true;
    } else if (i1 == j1) {
      if (i2 < j2) {
        // else neighbor's index is less than face's
        lt = true;
      } else if (i2 == j2) {

        int dx=0,dy=0,dz=0;
        
        if (rank_ >= 1) {
          ASSERT1 ("Face::operator <", "fx_ = %d",fx_,
                   (-1 <= fx_ && fx_ <= 1));
          ASSERT1 ("Face::operator <", "face.fx_ = %d",face.fx_,
                   (-1 <= face.fx_ && face.fx_ <= 1));
          dx = fx_-face.fx_+1;
        }
        
        if (rank_ >= 2) {
          ASSERT1 ("Face::operator <", "fy_ = %d",fy_,
                   (-1 <= fy_ && fy_ <= 1));
          ASSERT1 ("Face::operator <", "face.fy_ = %d",face.fy_,
                   (-1 <= face.fy_ && face.fy_ <= 1));
          dy = fy_-face.fy_+1;
        }
        
        if (rank_ >= 3) {
          ASSERT1 ("Face::operator <", "fz_ = %d",fz_,
                   (-1 <= fz_ && fz_ <= 1));
          ASSERT1 ("Face::operator <", "face.fz_ = %d",face.fz_,
                   (-1 <= face.fz_ && face.fz_ <= 1));
          dz = fz_-face.fz_+1;
        }
        
        // else subface is less than face's
        
        if (dx + 3*(dy + 3*(dz)) < 0) lt = true;
      }
    }
    return lt;
  }

  bool operator == (const Face & face)
  {
    const Index &i1 = this->index_block_;
    const Index &i2 = this->index_neighbor_;
    const Index &j1 = face.index_block_;
    const Index &j2 = face.index_neighbor_;

    const bool lx = (rank_ >= 1) ? (fx_ == face.fx_) : true;
    const bool ly = (rank_ >= 2) ? (fy_ == face.fy_) : true;
    const bool lz = (rank_ >= 3) ? (fz_ == face.fz_) : true;
    return ((((i1==j1) && (i2==j2)) ||
             ((i1==j2) && (i2==j1))) &&
            (lx && ly && lz));
  }

  bool operator != (const Face & face)
  { return ! (*this == face); }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    p | index_block_;
    p | index_neighbor_;
    p | ax_;
    p | ay_;
    p | az_;
    p | px_;
    p | py_;
    p | pz_;
    p | fx_;
    p | fy_;
    p | fz_;
    p | axis_;
    p | face_;
    p | rank_;
  }

  Index index_block () const
  { return index_block_; }

  Index index_neighbor () const
  { return index_neighbor_; }

  void get_subface (int *fx, int *fy=0, int *fz=0)
  {
    if (fx) (*fx) = (rank_ >= 1) ? fx_ : 0;
    if (fy) (*fy) = (rank_ >= 2) ? fy_ : 0;
    if (fz) (*fz) = (rank_ >= 3) ? fz_ : 0;
  }

  /// Return whether block and neighbor are adjacent for each axis
  void adjacency (bool *lx, bool *ly, bool *lz) const;

  /// Return the starting index along the face.  This is used for
  /// neighboring blocks in different mesh levels to determine the
  /// position of the fine block face in the coarse block face.  returns
  ///  ( [ 0 | bx/2 ], [0 | by/2 ], [0 | bz/2] ) 
  void get_offset (int *p_ix, int *p_iy, int *p_iz,
                   int bx, int by, int bz);

  /// For periodic boundaries with a single block, faces on the
  /// boundary are not unique (e.g. single Block could be any of 6
  /// facets) set_normal() is used to disambiguate the Face.  Error
  /// if not consistent
  void set_normal (int axis, int face)
  {
    ASSERT1 ("Face::set_normal", "axis = %d", axis,
             (0 <= axis && axis < rank_));
    ASSERT1 ("Face::set_normal", "face = %d", face,
             (face == -1 || face == +1));
    axis_ = axis;
    face_ = face;
  }

  /// Return axis associated with the Face (0, 1, 2)
  int axis() const
  { return axis_; }

  /// Return direction associated with the Face (-1 or 1)
  int face() const
  { return face_; }

  /// Return the dimensionality
  int rank() const
  { return rank_; }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  Index index_block_;
  Index index_neighbor_;
  /// Array size of array-of-octree mesh
  int ax_, ay_, az_;
  /// Whether respective axes are periodic
  bool px_, py_, pz_;
  /// Optional subface within a facet (used for MHD)
  int fx_,fy_,fz_;
  /// Axis associated with face normal.  Used to resolve ambiguity for
  /// periodic b.c. and single block along axis
  int axis_;
  /// Direction associated with face normal (-1 or 1).  Used to
  /// resolve ambiguity for periodic b.c. and single block along axis
  int face_;
  /// Dimensionality 1, 2, or 3
  int rank_;
};

#endif /* MESH_FACE_HPP */

