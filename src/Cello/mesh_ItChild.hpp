// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItChild.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    [\ref Mesh] Declaration of the ItChild class
///

#ifndef MESH_ITCHILD_HPP
#define MESH_ITCHILD_HPP

class ItChild {

  /// @class    ItChild
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  ItChild(int rank, const int * if3=0) throw()
    : ic3_(),
      rank_(rank)
  { reset();
    for (int i=0;     i<3; i++) 
      if3_[i] = (if3 && i<rank_) ? if3[i] : 0;
  }

  /// Destructor
  ~ItChild() throw()
  {
  }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    TRACEPUP;
    PUParray(p,ic3_,3);
    PUParray(p,if3_,3);
    p | rank_;
  }

  /// Go to next child, returning false when done
  bool next (int ic3[3]) throw()
  {

    do {
      increment_() ;
    } while (!valid_());

    ic3[0] = rank_ >= 1 ? ic3_[0] : 0;
    ic3[1] = rank_ >= 2 ? ic3_[1] : 0;
    ic3[2] = rank_ >= 3 ? ic3_[2] : 0;
    return (!is_reset());
  }

  /// Reset the Iterator to the beginning
  void reset() throw()
  {
    ic3_[0] = -1;
    ic3_[1] = 0;
    ic3_[2] = 0;
  };

  bool is_reset() const
  {
    return ic3_[0] == -1; 
  }

  /// Return the current value of the reduction operator
  void value(int ic3[3]) const throw()
  {
    ic3[0] = ic3_[0];
    ic3[0] = ic3_[0];
    ic3[0] = ic3_[0];
  };

private: // functions

  /// go to the next 

  void increment_()
  {
    if (is_reset()) {
      set_first_();
    } else {
      if (rank_ >= 1 && ic3_[0] < 1) {
	++ic3_[0];
      }	else {
	ic3_[0] = 0;
	if (rank_ >= 2 && ic3_[1] < 1) {
	  ++ic3_[1];
	} else {
	  ic3_[1] = 0;
	  if (rank_ >= 3 && ic3_[2] < 1) {
	    ++ic3_[2];
	  } else {
	    reset();
	  }
	}
      }
    }
  }

  /// Go to the first child
  void set_first_()
  {
    ic3_[0] = 0;
    ic3_[1] = 0;
    ic3_[2] = 0;
  }

  /// Whether the current face rank is valid
  bool valid_() const
  {
    if (is_reset()) return true;

    bool valid = true;
    if (rank_ >= 0 && if3_[0] == -1 && ic3_[0] != 0) valid = false;
    if (rank_ >= 0 && if3_[0] ==  1 && ic3_[0] != 1) valid = false;
    if (rank_ >= 1 && if3_[1] == -1 && ic3_[1] != 0) valid = false;
    if (rank_ >= 1 && if3_[1] ==  1 && ic3_[1] != 1) valid = false;
    if (rank_ >= 2 && if3_[2] == -1 && ic3_[2] != 0) valid = false;
    if (rank_ >= 2 && if3_[2] ==  1 && ic3_[2] != 1) valid = false;

    return valid;
  }

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Current child
  int ic3_[3];

  /// Adjacency face
  int if3_[3];

  /// simulation rank
  int rank_;

};

#endif /* MESH_ITCHILD_HPP */

