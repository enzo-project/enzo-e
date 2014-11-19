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
      field_ghosts_(field_ghosts),
      field_face_rank_(field_face_rank),
      field_list_()
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
    p | field_ghosts_;
    p | field_list_;
    p | field_face_rank_;
    p | name_;
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

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Ghost zone depth
  int field_ghosts_;

  /// Indicies of fields to include
  std::vector <int> field_list_;

  /// minimum face field rank to refresh (0 = corners, 1 = edges, etc.)
  int field_face_rank_;

  /// Name of this Refresh object
  std::string name_;

};

#endif /* PROBLEM_REFRESH_HPP */

