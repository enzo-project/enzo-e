// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Refresh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-11-04 22:24:46
/// @brief    [\ref Problem] Declaration of the Refresh class
///

#ifndef PROBLEM_REFRESH_HPP
#define PROBLEM_REFRESH_HPP

class Refresh : public PUP::able {

  /// @class    Refresh
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// empty constructor for charm++ pup()
  Refresh() throw() 
  : all_fields_(false),
    field_list_src_(),
    field_list_dst_(),
    all_particles_(false),
    particle_list_(),
    ghost_depth_(0),
    min_face_rank_(0),
    neighbor_type_(neighbor_leaf),
    accumulate_(false),
    sync_type_   (sync_unknown),
    sync_id_ (-1),
    active_(true),
    callback_(0) ,
    root_level_(0)
#ifdef NEW_REFRESH
    , id_refresh_(-1)
#endif
  {
  }

  /// Create an initialized Refresh object
  Refresh
  (int ghost_depth,
   int min_face_rank,
   int neighbor_type,
   int sync_type,
   int sync_id,
   bool active=true) throw() 
    : all_fields_(false),
      field_list_src_(),
      field_list_dst_(),
      all_particles_(false),
      particle_list_(),
      ghost_depth_(ghost_depth),
      min_face_rank_(min_face_rank),
      neighbor_type_(neighbor_type),
      accumulate_(false),
      sync_type_(sync_type),
      sync_id_(sync_id),
      active_(active),
      callback_(0),
      root_level_(0)
#ifdef NEW_REFRESH
    , id_refresh_(-1)
#endif      
  {
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(Refresh);

  /// CHARM++ migration constructor for PUP::able
  Refresh (CkMigrateMessage *m)
    : PUP::able(m),
    all_fields_(false),
    field_list_src_(),
    field_list_dst_(),
    all_particles_(false),
    particle_list_(),
    ghost_depth_(0),
    min_face_rank_(0),
    neighbor_type_(0),
    accumulate_(false),
    sync_type_(0),
    sync_id_ (-1),
    active_(true),
    callback_(0),
    root_level_(0)
#ifdef NEW_REFRESH    
    , id_refresh_(-1)
#endif    
  {
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    PUP::able::pup(p);

    p | all_fields_;
    p | field_list_src_;
    p | field_list_dst_;
    p | all_particles_;
    p | particle_list_;
    p | ghost_depth_;
    p | min_face_rank_;
    p | neighbor_type_;
    p | accumulate_;
    p | sync_type_;
    p | sync_id_;
    p | active_;
    p | callback_;
    p | root_level_;
#ifdef NEW_REFRESH
    p | id_refresh_;
#endif    
  }

  //--------------------------------------------------
  // FIELD METHODS
  //--------------------------------------------------

  /// Add a field id to the list of fields to refresh; don't add if it's already
  /// in the list
  void add_field(int id_field) {
    all_fields_ = false;
    if (id_field >= 0) {
      if (std::find (field_list_src_.begin(), field_list_src_.end(), id_field)
	  == field_list_src_.end()) {
	field_list_src_.push_back(id_field);
      }
      if (std::find (field_list_dst_.begin(), field_list_dst_.end(), id_field)
	  == field_list_dst_.end()) {
	field_list_dst_.push_back(id_field);
      }
    }
  }

  /// Add a named field to the list of fields to refresh
  void add_field(std::string field_name);

  /// Add a source and corresponding destination field to refresh;
  /// does not check if fields are already in the lists
  void add_field_src_dst(int id_field_src, int id_field_dst) {
    all_fields_ = false;
    if (id_field_src >= 0 && id_field_dst >= 0) {
      field_list_src_.push_back(id_field_src);
      field_list_dst_.push_back(id_field_dst);
    }
  }
  
  /// Add named fields to the list of fields to refresh
  void add_field_src_dst(std::string field_src, std::string field_dst);
  
  /// All fields are refreshed
  void add_all_fields() {
    all_fields_ = true;
  }

  /// Add specified fields
  void set_field_list (std::vector<int> field_list)
  {
    field_list_src_ = field_list;
    field_list_dst_ = field_list;
  }
  
  /// Return whether all fields are refreshed
  bool all_fields() const
  { return all_fields_; }

  /// Return whether any fields are refreshed
  bool any_fields() const
  { return (all_fields_ || (field_list_src_.size() > 0)); }

  /// Return the list of source fields participating in the Refresh operation
  std::vector<int> & field_list_src()
  { return field_list_src_; }

  /// Return the list of destination fields participating in the
  /// Refresh operation
  std::vector<int> & field_list_dst()
  { return field_list_dst_; }

  //--------------------------------------------------
  // PARTICLE METHODS
  //--------------------------------------------------

  /// Add the given particle type to the list
  void add_particle(int id_particle) {
    all_particles_ = false;
    particle_list_.push_back(id_particle);
  }

  /// All particles types are refreshed
  void add_all_particles() {
    all_particles_ = true;
  }

  /// Return whether all particles are refreshed
  bool all_particles() const
  { return all_particles_; }

  /// Return whether any particles are refreshed
  bool any_particles() const
  { return (all_particles_ || (particle_list_.size() > 0)); }
  
  /// Return the list of particles participating in the Refresh operation
  std::vector<int> & particle_list() {
    return particle_list_;
  }

  /// Add all data
  void add_all_data()
  {
    add_all_fields();
    add_all_particles();
  }

  /// Return whether there are any data to be refreshed
  bool any_data()
  { return (any_fields() || any_particles()); }
  
  /// Whether this particular Block is participating in the Refresh operation
  void set_active (bool active) 
  { active_ = active; }

  bool active () const 
  { return active_; }

  /// Callback function after refresh is completed
  int callback() const { return callback_; };

  /// Set the callback function for after the refresh operation is completed
  void set_callback(int callback) 
  { callback_ = callback; }

  /// Coarse level for neighbor_tree neighbor type
  int root_level() const { return root_level_; };

  /// Set the coarse level for  neighbor_tree neighbor type
  void set_root_level(int root_level) 
  { root_level_ = root_level; }

  /// Return the current minimum rank (dimension) of faces to refresh
  /// e.g. 0: everything, 1: omit corners, 2: omit corners and edges
  int min_face_rank() const 
  { return min_face_rank_; }

  /// Return the data field ghost depth
  int ghost_depth() const
  { return ghost_depth_; }

  /// Return the type of neighbors to refresh with: neighbor_leaf for
  /// neighboring leaf node (may be different mesh level) or
  /// neighbor_level for neighboring block in the same level (may be
  /// non-leaf)
  int neighbor_type() const 
  { return neighbor_type_; }

  /// Return whether to add neighbor face values to ghost zones or to
  /// copy them.  NOTE only accumulates if source field is different
  /// from destination field
  bool accumulate() const 
  {
    return accumulate_;
  }

  /// Set whether to add neighbor face values to ghost zones instead of
  /// copying them.
  void set_accumulate(bool accumulate)
  {
    accumulate_ = accumulate;
  }

  //----------------
  // Synchronization
  //----------------

  int sync_type() const 
  { return sync_type_; }

  int sync_load() const
  { return 3*sync_id_; }
  int sync_store() const
  { return 3*sync_id_+1; }
  int sync_exit() const
  { return 3*sync_id_+2; }

  void print() const 
  {
    CkPrintf ("Refresh %p\n",this);
    CkPrintf ("Refresh %p all_fields = %d\n",this,all_fields_);
    CkPrintf ("Refresh %p src fields:",this);
    for (size_t i=0; i<field_list_src_.size(); i++)
      CkPrintf (" %d",field_list_src_[i]);
    CkPrintf ("\n");
    CkPrintf ("Refresh %p dst fields:",this);
    for (size_t i=0; i<field_list_dst_.size(); i++)
      CkPrintf (" %d",field_list_dst_[i]);
    CkPrintf ("\n");
    CkPrintf ("Refresh %p all_particles = %d\n",this,all_particles_);
    CkPrintf ("Refresh %p particles:",this);
    for (size_t i=0; i<particle_list_.size(); i++)
      CkPrintf (" %d",particle_list_[i]);
    CkPrintf ("\n");
    CkPrintf ("Refresh %p ghost_depth = %d\n",this,ghost_depth_);
    CkPrintf ("Refresh %p min_face_rank: %d\n",this,min_face_rank_);
    CkPrintf ("Refresh %p neighbor_type: %d\n",this,neighbor_type_);
    CkPrintf ("Refresh %p accumulate: %d\n",this,accumulate_);
    CkPrintf ("Refresh %p sync_type: %d\n",this,sync_type_);
    CkPrintf ("Refresh %p sync_id: %d\n",this,sync_id_);
    CkPrintf ("Refresh %p active: %d\n",this,active_);
    CkPrintf ("Refresh %p callback: %d\n",this,callback_);
    CkPrintf ("Refresh %p root_level: %d\n",this,root_level_);
    fflush(stdout);
  }

  /// Return loop limits 0:3 for 4x4x4 particle data array indices
  /// for the given neighbor
  void get_particle_bin_limits
  (int rank,
   int refresh_type, 
   int if3[3], int ic3[3],
   int lower[3], int upper[3])
  {
    for (int axis=0; axis<rank; axis++) {
      if (if3[axis] == -1) { 
	lower[axis] = 0;
	upper[axis] = 1;
      } else if (if3[axis] == +1) { 
	lower[axis] = 3;
	upper[axis] = 4;
      } else {
	if (refresh_type == refresh_same) {
	  lower[axis] = 1;
	  upper[axis] = 3;
	} else if (refresh_type == refresh_fine) {
	  lower[axis] = ic3[axis] + 1;
	  upper[axis] = ic3[axis] + 2;
	} else if (refresh_type == refresh_coarse) {
	  lower[axis] = 1 - ic3[axis];
	  upper[axis] = 4 - ic3[axis];
	} else {
	  print();
	  ERROR1 ("Refresh::get_particle_bin_limits()",
		  "unknown refresh_type %d",
		  refresh_type);
	}
      }
    }
  }

#ifdef NEW_REFRESH
  /// Set the new refresh id in new_refresh_list_[]
  void set_id(int id_refresh)
  { id_refresh_ = id_refresh; }

  /// return the new refresh id in new_refresh_list_[]
  int id() const
  { return id_refresh_; }
#endif  
  
  //--------------------------------------------------

  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer);

  //--------------------------------------------------

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Whether to refresh all fields, ignoring field_list_*
  int all_fields_;
  
  /// Indicies of source fields; assumes all_fields_ == false, and
  /// size must be equal to field_list_dst_;
  std::vector <int> field_list_src_;

  /// Indicies of corresponding destination fields; assumes
  /// all_fields_ == false, and size must be equal to field_list_src_;
  std::vector <int> field_list_dst_;

  /// Whether to refresh all particle types, ignoring particle_list_
  int all_particles_;
  
  /// Whether to ignore particle_list_ and refresh all particle types
  /// Indicies of particles to include
  std::vector <int> particle_list_;

  /// Ghost zone depth
  int ghost_depth_;

  /// minimum face rank to refresh (2 include facets, 1 also edges, 0
  /// also corners)
  int min_face_rank_;

  /// Which subset of adjacent Blocks to refresh with
  int neighbor_type_;

  /// Whether to copy or add values
  int accumulate_;
  
  /// Synchronization type
  int sync_type_;

  /// Index for refresh synchronization counter
  int sync_id_;

  /// Whether the Refresh object is active for the block (replaces is_leaf())
  int active_;

  /// Callback after the refresh operation
  int callback_;

  /// Coarse level for neighbor_tree type
  int root_level_;

#ifdef NEW_REFRESH  
  /// ID in new_refresh_list_[]
  int id_refresh_;
#endif  
};

#endif /* PROBLEM_REFRESH_HPP */

