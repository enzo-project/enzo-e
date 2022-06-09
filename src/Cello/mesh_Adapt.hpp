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

  class LevelInfo {
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

  /// Constructor
  Adapt ()
    : valid_(false),
      rank_(0),
      min_level_(0),
      max_level_(0)
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
    PUParray(p,face_level_,3);

    p | valid_;
    p | rank_;
    PUParray (p,periodicity_,3);
    p | min_level_;
    p | max_level_;
    p | self_;
    p | neighbor_list_;

  }

  enum class LevelType { curr, next, last };

  //----------------------------------------------------------------------

  /// Set level for the block in the given face and (neighbor's) child
  void set_face_level (const int if3[3], LevelType level_type,
                       int level_min, int level_max, bool can_coarsen);

  /// Set level for the block in the given face and (neighbor's) child
  void set_face_level (const int if3[3], LevelType level_type,
                       int level_min)
  {
    face_level_[int(level_type)][IF3(if3)] = level_min;
  }

  /// Set level for the block in the given face and (neighbor's) child
  void set_face_level_last (const int ic3[3], const int if3[3],
                            int level_min, int level_max, bool can_coarsen);

  /// Set level for the block in the given face and (neighbor's) child
  void set_face_level_last (const int ic3[3], const int if3[3],
                            int level_min)
  {
    face_level_[int(LevelType::last)][ICF3(ic3,if3)] = level_min;
  }

  /// Return the current level_now for the block in the given face and
  /// (neighbor's) child
  int face_level (const int if3[3],LevelType level_type) const
  {
    ASSERT("Adapt::face_level",
           "face_level() called with LevelType::last",
           (level_type != LevelType::last));
    return face_level_[int(level_type)][IF3(if3)];
  }

  int face_level (int axis, int face, LevelType level_type) const
  {
    ASSERT("Adapt::face_level",
           "face_level() called with LevelType::last",
           (level_type != LevelType::last));
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
    ASSERT("Adapt::copy_face_level",
           "copy_face_level() called with LevelType::last",
           (level_type != LevelType::last));
    int n = face_level_[int(level_type)].size();
    for (int i=0; i<n; i++) face_level_[int(level_type)][i] = face_level[i];
  }

  inline void reset ()
  {
    reset_face_level(LevelType::curr);
    reset_face_level(LevelType::next);
    reset_face_level(LevelType::last);
  }

  //======================================================================

  /// Set rank. No range checking, rank must be 1 <= rank <= 3. This
  /// is only here so it can be called by test code outside of Cello,
  /// otherwise it would just use cello::rank()
  inline void set_rank (int rank)
  {
    rank_ = rank;
    set_valid(true);
  }

  inline void set_valid (bool valid)
  { valid_ = valid; }

  inline void set_periodicity (const int periodicity[3])
  {
    periodicity_[0] = periodicity[0];
    periodicity_[1] = periodicity[1];
    periodicity_[2] = periodicity[2];
  }
  /// Set the maximum allowable refinement level
  inline void set_max_level (int max_level)
  { max_level_ = max_level; }
  /// Set the minimum allowable refinement level
  inline void set_min_level (int min_level)
  { min_level_ = min_level; }


  /// Set the ith index, where i=0 is the block's own index
  inline void set_index(Index index)
  {
    self_.index_ = index;
    self_.level_now_ = index.level();
  }

  inline bool insert_neighbor (Index index)
  { return insert_neighbor (index,self_.index_.is_sibling(index)); }

  /// Insert the given neighbor into the list of neighbors. Return
  /// true if successful and false if neighbor already inserted
  bool insert_neighbor  (Index index, bool is_sibling);

  /// Delete the given neighbor from list of neighbors. Return true if
  /// successful and false if neighbor not found.
  bool delete_neighbor  (Index index);

  /// Reset self and level bounds for next adapt phase
  void reset_bounds();

  /// Replace the neighboring block with refined neighbors
  void refine_neighbor  (Index index);

  /// Replace the neighboring block with a coarsened neighbor. May
  /// delete any neighboring sibling blocks, and may be called
  /// separately for siblings
  void coarsen_neighbor  (Index index);

  /// Refine self, replacing blocks non-adjacent blocks with siblings
  void refine (const Adapt & adapt_parent, int ic3[3]);

  /// Coarsen self, replacing blocks non-adjacent blocks with siblings
  void coarsen (const Adapt & adapt_child);

  void initialize_self
  (Index index, int level_min, int level_now);

  /// Update a Block neighbor’s current level bounds, presumably from
  /// updated bounds received by the neighboring Block.
  void update_neighbor
  (Index index, int level_min, int level_max, bool can_coarsen);

  /// Update the Block’s own level bounds given the current list of
  /// neighbor level bounds. Returns true iff the values change, which
  /// can be used to determine whether or not to update its neighbors
  /// with new level bounds.
  bool update_bounds();

  /// Return whether the given Block (the default is the block itself)
  /// is converged; that is, whether its minimum and maximum level
  /// bounds are the same. This can be used to signal that the Block’s
  /// level is finally determined, and can thus call a global
  /// reduction.
  bool is_converged() const;

  /// Return whether all neighbors are converged. This is used to
  /// ensure call to contribute happens only after all expected
  /// neighbor messages have been received
  bool neighbors_converged() const;

  /// Return the current level bounds of the given Block (default is
  /// the block itself.)
  void get_level_bounds
  (int * level_min, int * level_max, bool * can_coarsen) const;

  /// Return the current level bounds for the specified neighbor
  void get_neighbor_level_bounds
  (Index index, int * level_min, int * level_max, bool * can_coarsen) const;

  /// Return the minimum level for the given block
  inline int level_min() const {return self_.level_min_; }
  inline int level_max() const {return self_.level_max_; }
  inline bool can_coarsen() const {return self_.can_coarsen_; }
  inline int num_neighbors() const { return neighbor_list_.size(); }
  inline bool is_sibling (int i) const { return neighbor_list_[i].is_sibling_; }
  inline Index index() const { return self_.index_; }
  inline Index index(int i) const { return neighbor_list_.at(i).index_; }

  /// Display Adapt attributes for debugging
  void print(std::string message, Block * block = nullptr) const;

  /// Write block's neighbors to file for debugging
  void write(std::string filename, Block * block, int cycle_start = 0);

  ///--------------------
  /// PACKING / UNPACKING
  ///--------------------
  
  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  char * load_data (char * buffer);

  /// Return whether the index is in the list of neighbors, and return
  /// its index if it is (or one past last index if not)
  bool is_neighbor (Index index, int * ip = 0) const;

private: // methods

  /// Return vector of neighbor indices (used to get copy before
  /// modifying list in loop)
  std::vector<Index> index_neighbors_() const;

  void copy_ (const Adapt & adapt);

  void delete_neighbor_ (int i);

  LevelInfo * neighbor_ (Index index)
  {
    const int n = neighbor_list_.size();
    for (int i=0; i<n; i++) {
      if (neighbor_list_[i].index_ == index)
        return &neighbor_list_[i];
    }
    return nullptr;
  }


private: // attributes

  /// List of neighbor indices (and self)
  // NOTE: change pup() function whenever attributes change
  std::vector<int> face_level_[3];

  /// Whether this Adapt class is valid; used for resetting existing
  /// Adapt
  bool valid_;
  /// Dimensionality of the problem
  int rank_;
  /// Periodicity
  int periodicity_[3];
  /// Minimum refinement level
  int min_level_;
  /// Maximum refinement level
  int max_level_;
  /// Level bound information for this block
  LevelInfo self_;
  /// Level bound information for neighboring blocks
  std::vector<LevelInfo> neighbor_list_;

};

#endif /* MESH_ADAPT_HPP */

