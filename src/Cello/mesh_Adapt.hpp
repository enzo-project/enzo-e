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
    : rank_(0)
  {
#ifdef OLD_ADAPT
    face_level_[0].resize(27);
    face_level_[1].resize(27);
    face_level_[2].resize(27*8);
    reset_face_level(LevelType::curr);
    reset_face_level(LevelType::next);
    reset_face_level(LevelType::last);
#endif
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
#ifdef OLD_ADAPT
    PUParray(p,face_level_,3);
#endif

    p | rank_;
    p | num_children_;
    p | level_now_;
    p | level_min_;
    p | level_max_;
    p | i_can_coarsen_;
    p | index_set_;
    p | index_map_;
    p | neighbor_list_;
    p | index_;

  }

  enum class LevelType { curr, next, last };

  //----------------------------------------------------------------------

#ifdef OLD_ADAPT
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

  /// Return the current level_now for the block in the given face and
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
#endif /* OLD_ADAPT */

  // OLD API
  //======================================================================
  // NEW API

  /// Set rank. No range checking, rank must be 1 <= rank <= 3
  inline void set_rank (int rank)
  { rank_ = rank; }

  /// Set the ith index, where i=0 is the block's own index
  inline void set_index(Index index)
  { index_ = index; }

  inline bool insert_neighbor (Index index)
  { return insert_neighbor (index,index_.is_sibling(index)); }

  /// Insert the given neighbor into the list of neighbors. Return
  /// true if successful and false if neighbor already inserted
  bool insert_neighbor  (Index index, bool is_sibling);

  /// Delete the given neighbor from list of neighbors. Return true if
  /// successful and false if neighbor not found.
  bool delete_neighbor  (Index index);

  /// Replace the neighboring block with refined neighbors
  bool refine_neighbor  (Index index);

  /// Replace the neighboring block with a coarsened neighbor. May
  /// delete any neighboring sibling blocks, and may be called
  /// separately for siblings
  bool coarsen_neighbor  (Index index);

  /// Refine self, replacing blocks non-adjacent blocks with siblings
  void refine (const Adapt & adapt_parent, int ic3[3]);

  /// Coarsen self, replacing blocks non-adjacent blocks with siblings
  void coarsen (const Adapt & adapt_child);

  void initialize_self
  (Index index, int level_min, int level_now, int level_max);

  /// Update a Block neighbor’s current level bounds, presumably from
  /// updated bounds received by the neighboring Block.
  void update_neighbor
  (Index index, int level_min, int level_max, bool can_coarsen);

  /// Update the Block’s own level bounds given the current list of
  /// neighbor level bounds. Returns true iff the values change, which
  /// can be used to determine whether or not to update its neighbors
  /// with new level bounds.
  bool update_bounds ();

  /// Return whether the given Block (the default is the block itself)
  /// is “committed”; that is, whether its minimum and maximum level
  /// bounds are the same. This can be used to signal that the Block’s
  /// level is finally determined, and can thus call a global
  /// reduction.
  bool is_committed() const;

  /// Return the current level bounds of the given Block (default is
  /// the block itself.)
  void get_level_bounds
  (int * level_min, int * level_max, bool * can_coarsen) const;

  /// Return the element of the vector storing the neighboring Block
  /// level bounds, where 0 is reserved for the Block itself.
  //  inline int index (Index index) const { return index_map_.at(index); };

  //   inline Index index (int i) const { return index_list_.at(i); }

  /// Return the minimum level for the given block
  inline int level_min() const {return level_min_; }
  inline int level_max() const {return level_max_; }
  inline bool can_coarsen() const {return i_can_coarsen_; }
  inline int num_neighbors() const { return neighbor_list_.size(); }
  inline bool is_sibling (int i) const { return neighbor_list_[i].is_sibling_; }
  inline Index index() const { return index_; }
  inline Index index(int i) const { return neighbor_list_.at(i).index_; }
  void print(std::string message) const;

private: // methods

  /// Return whether the index is in the list of neighbors, and return
  /// its index if it is (or one past last indicex if not)
  bool is_neighbor_ (Index index, int * ip = 0) const;

  void copy_ (const Adapt & adapt);

  void delete_neighbor_ (int i);
  
private: // attributes

#ifdef OLD_ADAPT
  // NOTE: change pup() function whenever attributes change
  std::vector<int> face_level_[3];
#endif

  /// Dimensionality of the problem
    int rank_;
  
    /// Number of children per block; set using cello::num_children() but
    /// may be updated for testing via set_rank
    int num_children_;

    int level_now_;
    int level_min_;
    int level_max_;

    /// Current (starting) level for neighbors
    //  std::vector<int> level_now_;

    /// Level bounds for neighbors
    //  std::vector<int> level_min_;
    //  std::vector<int> level_max_;

    /// Whether the neighbor block is a sibling (shares same coarse
    /// parent)
    //  std::vector<bool> is_sibling_;

    /// Whether the sibling block could coarsen if all of its siblings
    /// can coarsen.
    //  std::vector<bool> can_coarsen_;

    bool i_can_coarsen_;

    /// Mapping of neighbor (and self) index to level bound arrays
    std::set<Index> index_set_;
    std::map<Index,int> index_map_;

  /// List of neighbor indices (and self)
  class NeighborInfo {
  public:
    void pup (PUP::er &p)
    {
      p | index_;
      p | level_now_;
      p | level_min_;
      p | level_max_;
      p | is_sibling_;
      p | can_coarsen_;
    }
    Index index_;
    int level_now_;
    int level_min_;
    int level_max_;
    bool is_sibling_;
    bool can_coarsen_;
  };
    
  std::vector<NeighborInfo> neighbor_list_;

  //  std::vector<Index>  index_list_;
  Index index_;
};

#endif /* MESH_ADAPT_HPP */

