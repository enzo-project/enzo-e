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
  : field_list_(),
    particle_list_(),
    ghost_depth_(0),
    min_face_rank_(0),
    neighbor_type_(neighbor_unknown),
    sync_type_   (sync_unknown),
    sync_load_(),
    sync_store_(),
    sync_id_ (0),
    active_(true),
    callback_(0) 
  {
  }

  /// empty constructor for charm++ pup()
  Refresh
  (int ghost_depth,
   int min_face_rank,
   int neighbor_type,
   int sync_type,
   int sync_id=0,
   bool active=true) throw() 
  : field_list_(),
    particle_list_(),
    ghost_depth_(ghost_depth),
    min_face_rank_(min_face_rank),
    neighbor_type_(neighbor_type),
    sync_type_(sync_type),
    sync_load_(),
    sync_store_(),
    sync_id_(sync_id),
    active_(active),
    callback_(0) 
  {
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(Refresh);

  /// CHARM++ migration constructor for PUP::able
  Refresh (CkMigrateMessage *m)
    : PUP::able(m),
      field_list_(),
      particle_list_(),
      ghost_depth_(0),
      min_face_rank_(0),
      neighbor_type_(0),
      sync_type_(0),
      sync_load_(),
      sync_store_(),
      sync_id_ (0),
      active_(false),
      callback_(0)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    PUP::able::pup(p);

    p | field_list_;
    p | particle_list_;
    p | ghost_depth_;
    p | min_face_rank_;
    p | neighbor_type_;
    p | sync_type_;
    p | sync_load_;
    p | sync_store_;
    p | sync_id_;
    p | active_;
    p | callback_;
  }

  //--------------------------------------------------
  // FIELD METHODS
  //--------------------------------------------------

  /// Add the given field to the list
  void add_field(int id_field) {
    field_list_.push_back(id_field);
  }

  /// All fields are refreshed
  void add_all_fields(int num_fields) {
    field_list_.clear();
    for (int i=0; i<num_fields; i++) {
      field_list_.push_back(i);
    }
  }

  /// Return the list of fields participating in the Refresh operation
  std::vector<int> & field_list() {
    return field_list_;
  }

  //--------------------------------------------------
  // PARTICLE METHODS
  //--------------------------------------------------

  /// Add the given particle type to the list
  void add_particle(int id_particle) {
    particle_list_.push_back(id_particle);
  }

  /// All particles types are refreshed
  void add_all_particles(int num_particles) {
    particle_list_.clear();
    for (int i=0; i<num_particles; i++) {
      particle_list_.push_back(i);
    }
  }

  /// Return the list of particles participating in the Refresh operation
  std::vector<int> & particle_list() {
    return particle_list_;
  }
  /// Whether this Block is participating in the Refresh operation

  void set_active (bool active) 
  { active_ = active; }

  bool active () const 
  { return active_; }

  /// Callback function after refresh is completed
  int callback() const { return callback_; };

  void set_callback(int callback) 
  { callback_ = callback; }

  /// Return the current minimum rank (dimension) of faces to
  /// refresh e.g. 0: everything, 1: omit corners, 2: omit corners and
  /// edges

  int min_face_rank() const 
  { return min_face_rank_; }

  /// Return the data field ghost depth

  int ghost_depth() const
  { return ghost_depth_; }

  int neighbor_type() const 
  { return neighbor_type_; }

  int sync_type() const 
  { return sync_type_; }

  Sync & sync_load() 
  {  return sync_load_; }

  Sync & sync_store() 
  {  return sync_store_; }

  int sync_id() 
  {  return sync_id_; }


  void print() const 
  {
    CkPrintf ("Refresh %p\n",this);
    CkPrintf ("Refresh fields:");
    for (size_t i=0; i<field_list_.size(); i++)
      CkPrintf (" %d",field_list_[i]);
    CkPrintf ("\n");
    CkPrintf ("Refresh particles:");
    for (size_t i=0; i<particle_list_.size(); i++)
      CkPrintf (" %d",particle_list_[i]);
    CkPrintf ("\n");
    CkPrintf ("Refresh ghost_depth = %d\n",ghost_depth_);
    CkPrintf ("Refresh min_face_rank: %d\n",min_face_rank_);
    CkPrintf ("Refresh neighbor_type: %d\n",neighbor_type_);
    CkPrintf ("Refresh sync_type: %d\n",sync_type_);
    CkPrintf ("Refresh sync_load: %d/%d\n", 
	    sync_load_.value(),
	    sync_load_.stop());
    CkPrintf ("Refresh sync_store: %d/%d\n",
	    sync_store_.value(),
	    sync_store_.stop());
    CkPrintf ("Refresh sync_id: %d\n",sync_id_);
    CkPrintf ("Refresh active: %d\n",active_);
    CkPrintf ("Refresh callback: %d\n",callback_);
    fflush(stdout);
  }

  /// Return loop limits 0:3 for 4x4x4 particle data array indices
  /// for the given neighbor
  void index_limits
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
	  ERROR1 ("Refresh::loop_limits()",
		  "unknown refresh_type %d",
		  refresh_type);
	}
      }
    }
  }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Indicies of fields to include
  std::vector <int> field_list_;

  /// Indicies of particles to include
  std::vector <int> particle_list_;

  /// Ghost zone depth
  int ghost_depth_;

  /// minimum face rank to refresh (0 = corners, 1 = edges, etc.)
  int min_face_rank_;

  /// Which subset of adjacent Blocks to refresh with
  int neighbor_type_;

  /// Synchronization type
  int sync_type_;

  /// Counter for synchronization before loading data
  Sync sync_load_;

  /// Counter for synchronization after storing data
  Sync sync_store_;

  /// Index for refresh synchronization counter
  int sync_id_;

  /// Whether the Refresh object is active for the block (replaces is_leaf())
  bool active_;

  /// Callback after the refresh operation
  int callback_;

};

#endif /* PROBLEM_REFRESH_HPP */

