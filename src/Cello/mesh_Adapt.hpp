// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Adapt.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-10-20
/// @brief    [\ref Mesh] Declaration of the Adapt class

// #define CHECK_ADAPT

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
    void clear() {
      level_now_ = 0;
      level_min_ = 0;
      level_max_ = 0;
      is_sibling_ = false;
      can_coarsen_ = false;
    }
  };

  /// Constructor
  Adapt ()
    : face_level_curr_(),
      face_level_next_(),
      face_level_last_(),
      face_level_curr_count_(),
      face_level_next_count_(),
      face_level_last_count_(),
      valid_(false),
      rank_(0),
      periodicity_(),
      min_level_(0),
      max_level_(0),
      neighbor_list_(),
      count_(0)
  {
    reset_face_level_curr();
    reset_face_level_next();
    reset_face_level_last();
    periodicity_[0] = 0;
    periodicity_[1] = 0;
    periodicity_[2] = 0;
    self_.clear();
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    p | face_level_curr_;
    p | face_level_next_;
    p | face_level_last_;
    p | face_level_curr_count_;
    p | face_level_next_count_;
    p | face_level_last_count_;
    p | valid_;
    p | rank_;
    PUParray (p,periodicity_,3);
    p | min_level_;
    p | max_level_;
    p | self_;
    p | neighbor_list_;
    p | count_;

  }

  //----------------------------------------------------------------------

  /// Set current level for the block in the given face if counter is
  /// greater than saved
  void set_face_level_curr (const int if3[3], int level, int count);

  /// Set next level for the block in the given face if counter is
  /// greater than saved
  void set_face_level_next (const int if3[3], int level, int count);

  /// Set last level for the block in the given face and neighbor's
  /// child if counter is greater than last
  void set_face_level_last
  (const int ic3[3], const int if3[3], int level, int count);

  /// Return the current level_now for the block in the given face and
  /// (neighbor's) child
  int face_level_curr (const int if3[3]) const
  { return face_level_curr_[IF3(if3)];  }

  int face_level_next (const int if3[3]) const
  { return face_level_next_[IF3(if3)];  }

  int face_level_curr (int axis, int face) const
  {
    int if3[3];
    cello::af_to_xyz(axis,face,if3);
    return face_level_curr(if3);
  }

  int face_level_next (int axis, int face) const
  {
    int if3[3];
    cello::af_to_xyz(axis,face,if3);
    return face_level_next(if3);
  }

  int face_level_last (const int ic3[3], const int if3[3]) const
  {
    return face_level_last_[IC3(ic3)][IF3(if3)];
  }

  inline void reset_face_level_curr()
  {
    clear_curr_face_count();
    if (face_level_curr_.size() != 27)
        face_level_curr_.resize(27);
    std::fill(face_level_curr_.begin(),
              face_level_curr_.end(), 0 );
  }

  inline void reset_face_level_next()
  {
    clear_next_face_count();
    if (face_level_next_.size() != 27)
        face_level_next_.resize(27);
    std::fill(face_level_next_.begin(),
              face_level_next_.end(), 0 );
  }

  inline void reset_face_level_last()
  {
    clear_last_face_count();
    const int nc = cello::num_children();
    if (face_level_last_.size() != nc)
        face_level_last_.resize(nc);
    for (size_t i=0; i<nc; i++) {
      if (face_level_last_[i].size() != 27)
          face_level_last_[i].resize(27);
      std::fill(face_level_last_[i].begin(),
                face_level_last_[i].end(), -1 );
    }
  }

  /// Clear all face level counts
  void clear_face_counts()
  {
    clear_curr_face_count();
    clear_next_face_count();
    clear_last_face_count();
  }

  /// Clear current face level counts
  void clear_curr_face_count()
  {
    auto & curr = face_level_curr_count_;
    if (curr.size() != 27) curr.resize(27);
    std::fill(curr.begin(), curr.end(), -1 );
  }

  /// clear next face level counts
  void clear_next_face_count()
  {
    auto & next = face_level_next_count_;
    if (next.size() != 27) next.resize(27);
    std::fill(next.begin(), next.end(), -1 );
  }

  /// clear last face level counts
  void clear_last_face_count()
  {
    auto & last = face_level_last_count_;
    const int nc = cello::num_children();
    if (last.size() != nc) last.resize(nc);
    for (size_t i=0; i<nc; i++) {
      if (last[i].size() != 27) last[i].resize(27);
      std::fill(last[i].begin(), last[i].end(), -1 );
    }
  }

  /// Return vector of current face levels
  std::vector<int> & face_level_curr()
  { return face_level_curr_; }

  /// Return vector of current face level counts
  std::vector<int> & face_level_curr_count()
  { return face_level_curr_count_; }

  void update_curr_from_next ()
  { face_level_curr_ = face_level_next_; }

  void update_next_from_curr ()
  { face_level_next_ = face_level_curr_; }

  void copy_face_level_curr(int * face_level)
  {
    for (int i=0; i<27; i++)
      face_level_curr_[i] = face_level[i];
  }

  inline void reset ()
  {
    reset_face_level_curr();
    reset_face_level_next();
    reset_face_level_last();
  }

  //======================================================================

  /// Set rank. No range checking, rank must be 1 <= rank <= 3
  inline void set_rank (int rank)
  {
    rank_ = rank;
    set_valid(true);
  }

  inline void set_valid (bool valid)
  { valid_ = valid; }

  inline bool is_valid () const
  { return valid_; }

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

  inline void insert_neighbor (Index index)
  { insert_neighbor (index,self_.index_.is_sibling(index)); }

  /// Insert the given neighbor into the list of neighbors. Return
  /// true if successful and false if neighbor already inserted
  void insert_neighbor  (Index index, bool is_sibling);

  /// Delete the given neighbor from list of neighbors. Return true if
  /// successful and false if neighbor not found.
  bool delete_neighbor  (Index index);

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

  /// Update a Block neighbor’s current level bounds, presumably from
  /// updated bounds received by the neighboring Block.
  void update_neighbor
  (Index index, int level_min, int level_max, bool can_coarsen);

  /// Initialize level bounds for self and neighbors
  void initialize_bounds(int level_next, int level_curr);

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

  /// Return vector of neighbor indices (used to get copy before
  /// modifying list in loop)
  std::vector<Index> index_neighbors() const;

  /// Return the minimum level for the given block
  inline int level_min() const {return self_.level_min_; }
  inline int level_max() const {return self_.level_max_; }
  inline int level_now() const {return self_.level_now_; }
  inline bool can_coarsen() const {return self_.can_coarsen_; }
  inline int num_neighbors() const { return neighbor_list_.size(); }
  inline bool is_sibling (int i) const { return neighbor_list_[i].is_sibling_; }
  inline Index index() const { return self_.index_; }
  inline Index index(int i) const { return neighbor_list_.at(i).index_; }

  /// Display Adapt attributes for debugging
  void print(std::string message, const Block * block = nullptr, FILE * fp=nullptr) const;

  /// Write block's neighbors to file for debugging
  void write(std::string filename, const Block * block, int cycle_start = 0) const;

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

  int count() const { return count_; }
  void inc_count () { ++count_; }
  void clear_count() { count_ = 0; }

private: // methods

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

  // NOTE: change pup() function whenever attributes change

  /// List of neighbor levels
  std::vector<int> face_level_curr_;
  std::vector<int> face_level_next_;
  std::vector < std::vector<int> > face_level_last_;

  /// List of counts associated with last assignment to neighbor levels
  /// Used to maintain assignment order even when Charm++ messages don't
  /// arrive in the same order sent
  std::vector<int> face_level_curr_count_;
  std::vector<int> face_level_next_count_;
  std::vector < std::vector<int> > face_level_last_count_;

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
  /// Local count
  int count_;
};

#endif /* MESH_ADAPT_HPP */

