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
  Refresh() throw() {}

  /// Create a Refresh object
  Refresh(std::string name,
	  int field_ghosts,
	  int field_face_rank) throw()
    : name_(name),
      field_list_(),
      field_ghosts_(field_ghosts),
      field_face_rank_(field_face_rank)
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

    p | name_;
    p | field_list_;
    p | field_ghosts_;
    p | field_face_rank_;
  }

  /// Add the given field to the list
  void insert_field(int id_field) {
    field_list_.push_back(id_field);
  }

  /// Set all fields
  void all_fields(int num_fields) {
    field_list_.clear();
    for (int i=0; i<num_fields; i++) {
      field_list_.push_back(i);
    }
  }

  std::vector<int> & field_list() {
    return field_list_;
  }

  /// Return the current minimum rank (dimension) of faces to refresh
  /// e.g. 0: everything, 1: omit corners, 2: omit corners and edges
  int field_face_rank() const 
  { return field_face_rank_; }

  /// Return the ghost zone depth
  int field_ghosts() const
  { return field_ghosts_; }

  /// Return the ith field index, or false if i is out of range
  bool get_field_index (size_t i, int * index_field)
  {
    const bool in_range = (i < field_list_.size());
    if (in_range) (*index_field) = field_list_[i];
    return in_range;
  }
  void print () const {
    printf ("%s:%d\n",__FILE__,__LINE__);
    printf ("name_ %s\n",name_.c_str());
    printf ("field_list:");
    for (size_t i=0; i<field_list_.size(); i++) {
      printf (" %d",field_list_[i]);
    }
    printf ("\n");
    printf ("field_ghosts: %d\n",field_ghosts_);
    printf ("field_face_rank: %d\n",field_face_rank_);

  }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Name of this Refresh object
  std::string name_;

  /// Indicies of fields to include
  std::vector <int> field_list_;

  /// Ghost zone depth
  int field_ghosts_;

  /// minimum face field rank to refresh (0 = corners, 1 = edges, etc.)
  int field_face_rank_;

};

#endif /* PROBLEM_REFRESH_HPP */

