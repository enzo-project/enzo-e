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
  Refresh() throw() { }

  /// Create a Refresh object
  Refresh(int ghost_depth,
	  int min_face_rank) throw()
    : field_list_(),
      ghost_depth_(ghost_depth),
      min_face_rank_(min_face_rank),
      sync_(),
      is_active_(true)

  {  }

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
    p | sync_type_;
    p | sync_;
    p | is_active_;
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

  /// Return whether specific field is included in refresh list
  bool field (int id_field) const {
    for (size_t i=0; i<field_list_.size(); i++) {
      if (field_list_.at(i) == id_field) return true;
    }
    return false;
  }
  std::vector<int> & field_list() {
    return field_list_;
  }

  void set_active (bool active) 
  { is_active_ = active; }
    
  /// Return the current minimum rank (dimension) of faces to refresh
  /// e.g. 0: everything, 1: omit corners, 2: omit corners and edges
  int min_face_rank() const 
  { return min_face_rank_; }

  /// Set ghost depth
  void set_ghost_depth(int ghost_depth)
  { ghost_depth_ = ghost_depth; }

  /// Return the ghost zone depth
  int ghost_depth() const
  { return ghost_depth_; }

  /// Return the ith field index, or false if i is out of range
  bool get_field_index (size_t i, int * index_field)
  {
    const bool in_range = (i < field_list_.size());
    if (in_range) (*index_field) = field_list_[i];
    return in_range;
  }

  void print () const {
    printf ("%s:%d\n",__FILE__,__LINE__);
    printf ("field_list:");
    for (size_t i=0; i<field_list_.size(); i++) {
      printf (" %d",field_list_[i]);
    }
    printf ("\n");
    printf ("ghost_depth: %d\n",ghost_depth_);
    printf ("min_face_rank: %d\n",min_face_rank_);

  }

  Sync & sync() 
  {  return sync_; }

  void set_sync_type(std::string sync_type) 
  { sync_type_ = sync_type; }

  std::string sync_type() const 
  { return sync_type_; }

  int callback() const { return callback_; };

  void set_callback(int callback) 
  { callback_ = callback; }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Indicies of fields to include
  std::vector <int> field_list_;

  /// Ghost zone depth
  int ghost_depth_;

  /// minimum face field rank to refresh (0 = corners, 1 = edges, etc.)
  int min_face_rank_;

  /// Synchronization type
  std::string sync_type_;

  /// Counter for synchronization
  Sync sync_;

  /// Whether the Refresh object is active for the block
  bool is_active_;

  /// Callback after the refresh operation
  int callback_;
};

#endif /* PROBLEM_REFRESH_HPP */

