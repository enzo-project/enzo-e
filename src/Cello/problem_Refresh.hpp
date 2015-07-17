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
    ghost_depth_(0),
    min_face_rank_(0),
    neighbor_type_(neighbor_unknown),
    sync_type_   (sync_unknown),
    sync_load_(),
    sync_store_(),
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
   bool active=true) throw() 
  : field_list_(),
    ghost_depth_(ghost_depth),
    min_face_rank_(min_face_rank),
    neighbor_type_(neighbor_type),
    sync_type_(sync_type),
    sync_load_(),
    sync_store_(),
    active_(active),
    callback_(0) 
  {
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(Refresh);

  /// CHARM++ migration constructor for PUP::able
  Refresh (CkMigrateMessage *m) : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    PUP::able::pup(p);

    p | field_list_;
    p | ghost_depth_;
    p | min_face_rank_;
    p | neighbor_type_;
    p | sync_type_;
    p | sync_load_;
    p | sync_store_;
    p | active_;
    p | callback_;
  }

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


  void print() const 
  {
    printf ("Refresh %p\n",this);
    printf ("Refresh fields:");
    for (size_t i=0; i<field_list_.size(); i++)
      printf (" %d",field_list_[i]);
    printf ("\n");
    printf ("Refresh ghost_depth = %d\n",ghost_depth_);
    printf ("Refresh min_face_rank: %d\n",min_face_rank_);
    printf ("Refresh neighbor_type: %d\n",neighbor_type_);
    printf ("Refresh sync_type: %d\n",sync_type_);
    printf ("Refresh sync_load: %d/%d\n", 
	    sync_load_.value(),
	    sync_load_.stop());
    printf ("Refresh sync_store: %d/%d\n",
	    sync_store_.value(),
	    sync_store_.stop());
    printf ("Refresh active: %d\n",active_);
    printf ("Refresh callback: %d\n",callback_);
  }
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Indicies of fields to include
  std::vector <int> field_list_;

  /// Ghost zone depth
  int ghost_depth_;

  /// minimum face field rank to refresh (0 = corners, 1 = edges, etc.)
  int min_face_rank_;

  /// Which subset of adjacent Blocks to refresh with
  int neighbor_type_;

  /// Synchronization type
  int sync_type_;

  /// Counter for synchronization before loading data
  Sync sync_load_;

  /// Counter for synchronization after storing data
  Sync sync_store_;

  /// Whether the Refresh object is active for the block (replaces is_leaf())
  bool active_;

  /// Callback after the refresh operation
  int callback_;

};

#endif /* PROBLEM_REFRESH_HPP */

