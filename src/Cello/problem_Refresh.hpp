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

  /// CHARM++ PUP::able declaration
  PUPable_decl(Refresh);

  /// CHARM++ migration constructor for PUP::able
  Refresh (CkMigrateMessage *m) : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { PUP::able::pup(p);
    TRACEPUP;
    p | field_set_;
    p | ghost_depth_;
  }

  /// Add an item to a group
  void add_field(int index_field)
    throw(std::out_of_range)
  {
    field_set_.insert(index_field);
  }

  /// Return list of Field's to Refresh
  const bool field(int index_field) const
  { return (field_set_.find(index_field) != field_set_.end()); }

  /// Set the ghost zone depth
  void set_ghost_depth(int ghost_depth)
  { ghost_depth_ = ghost_depth; }

  /// Return ghost zone depth
  int ghost_depth() const
  { return ghost_depth_; }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  std::set<int> field_set_;

  int ghost_depth_;

};

#endif /* PROBLEM_REFRESH_HPP */

