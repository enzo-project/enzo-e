// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Adapt.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-10-20
/// @brief    [\ref Mesh] Declaration of the Adapt class

#ifndef MESH_ADAPT_HPP
#define MESH_ADAPT_HPP

class Adapt {

  /// @class    Adapt
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh]

public: // interface

  /// Constructor
  Adapt ()
    : num_children_(cello::num_children())
  {
    face_level_[0].resize(27);
    face_level_[1].resize(27);
    face_level_[2].resize(27*8);
    reset_face_level(LevelType::curr);
    reset_face_level(LevelType::next);
    reset_face_level(LevelType::last);
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    p | face_level_[0];
    p | face_level_[1];
    p | face_level_[2];

    p | num_children_;
    p | index_map_;
    p | level_curr_;
    p | level_want_;
    p | level_min_;
    p | level_max_;
    p | is_sibling_;
    p | can_coarsen_;

  }

  enum class LevelType { curr, next, last };

  //----------------------------------------------------------------------

  /// Set level for the block in the given face and (neighbor's) child
  void set_face_level (const int if3[3], LevelType level_type,int level)
  {
    ASSERT("Adapt::set_face_level","set_face_level() called with LevelType::last",(level_type != LevelType::last));
    face_level_[int(level_type)][IF3(if3)] = level;
  }

  /// Set level for the block in the given face and (neighbor's) child
  void set_face_level_last (const int ic3[3], const int if3[3], int level)
  {
    face_level_[int(LevelType::last)][ICF3(ic3,if3)] = level;
  }

  /// Return the current level_curr for the block in the given face and
  /// (neighbor's) child
  int face_level (const int if3[3],LevelType level_type) const
  {
    ASSERT("Adapt::face_level","face_level() called with LevelType::last",(level_type != LevelType::last));
    return face_level_[int(level_type)][IF3(if3)];
  }

  int face_level (int axis, int face, LevelType level_type) const
  {
    ASSERT("Adapt::face_level","face_level() called with LevelType::last",(level_type != LevelType::last));
    int if3[3];
    cello::af_to_xyz(axis,face,if3);
    return face_level(if3,level_type);
  }

  int face_level_last (const int ic3[3], const int if3[3]) const
  {
    return face_level_[int(LevelType::last)][ICF3(ic3,if3)];
  }

  inline void reset_face_level (LevelType level_type)
  {
    int value = (level_type == LevelType::last) ? -1 : 0;
      std::fill(face_level_[int(level_type)].begin(),
                face_level_[int(level_type)].end(),value);
  }

  size_t size_face_level(LevelType level_type)
  { return face_level_[int(level_type)].size(); }

  std::vector<int> & vector_face_level(LevelType level_type)
  { return face_level_[int(level_type)]; }

  void update_curr_from_next ()
  {
    face_level_[int(LevelType::curr)] = face_level_[int(LevelType::next)];
  }
  void update_next_from_curr ()
  {
    face_level_[int(LevelType::next)] = face_level_[int(LevelType::curr)];
  }

  void copy_face_level(LevelType level_type, int * face_level)
  {
    ASSERT("Adapt::copy_face_level","copy_face_level() called with LevelType::last",(level_type != LevelType::last));
    int n = face_level_[int(level_type)].size();
    for (int i=0; i<n; i++) face_level_[int(level_type)][i] = face_level[i];
  }


  inline void reset ()
  {
    reset_face_level(LevelType::curr);
    reset_face_level(LevelType::next);
    reset_face_level(LevelType::last);
  }

  // OLD API
  //======================================================================
  // NEW API

  /// Set rank. No range checking, rank must be 1 <= rank <= 3
  void set_rank (int rank)
  {
    num_children_ = (rank == 1) ? 2 : (rank == 2) ? 4 : 8;
  }

  /// (Re)allocate storage for level bounds assuming num_neighbors
  /// adjacent Blocks. Includes extra element for the Block’s own
  /// level bounds
  void allocate_level_bounds (int num_neighbors);

  /// Initialize the level bounds for the current Block. Index i=0 is
  /// reserved for the Block itself; other indices are arbitrary, but
  /// typically ordered according to the ItNeighbor iterator.
  void set_initial_level_bounds
  (int i, Index index, int level_curr, int level_want, bool is_sibling);

  /// Update a Block neighbor’s current level bounds, presumably from
  /// updated bounds received by the neighboring Block.
  void update_level_bounds
  (Index index, int level_min, int level_max, bool can_coarsen);

  /// Update the Block’s own level bounds given the current list of
  /// neighbor level bounds. Returns true iff the values change, which
  /// can be used to determine whether or not to update its neighbors
  /// with new level bounds.
  bool evaluate_level_bounds ();

  /// Return whether the given Block (the default is the block itself)
  /// is “committed”; that is, whether its minimum and maximum level
  /// bounds are the same. This can be used to signal that the Block’s
  /// level is finally determined, and can thus call a global
  /// reduction.
  bool is_committed(int i = 0) const;

  /// Return the current level bounds of the given Block (default is
  /// the block itself.)
  void get_level_bounds
  (int * level_min, int * level_max, bool * can_coarsen, int i = 0) const;

  /// Return the element of the vector storing the neighboring Block
  /// level bounds, where 0 is reserved for the Block itself.
  int index (Index index) const { return index_map_.at(index); };

  /// Return the minimum level for the given block
  int level_curr(int i=0) const {return level_curr_.at(i); }
  int level_want(int i=0) const {return level_want_.at(i); }
  int level_min(int i=0) const {return level_min_.at(i); }
  int level_max(int i=0) const {return level_max_.at(i); }
  bool can_coarsen(int i=0) const {return can_coarsen_.at(i); }
  int num_neighbors() const { return level_min_.size(); }
  bool is_sibling (int i) const { return is_sibling_[i]; }
private: // attributes

  // NOTE: change pup() function whenever attributes change
  std::vector<int> face_level_[3];

  /// Number of children per block; set using cello::num_children() but
  /// may be updated for testing via set_rank
  int num_children_;
  
  /// Level bounds for self [0] and neighbors [1...n]
  std::vector<int> level_curr_;
  std::vector<int> level_want_;
  std::vector<int> level_min_;
  std::vector<int> level_max_;

  /// Whether the neighbor block is a sibling (shares same coarse
  /// parent)
  std::vector<bool> is_sibling_;

  /// Whether the sibling block could coarsen if all of its siblings
  /// can coarsen.
  std::vector<bool> can_coarsen_;

  /// Mapping of neighbor (and self) index to level bound arrays
  std::map<Index,int> index_map_;
};

#endif /* MESH_ADAPT_HPP */

