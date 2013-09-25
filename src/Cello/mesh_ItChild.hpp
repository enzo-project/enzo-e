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
  ItChild(int rank_simulation) throw()
    : ic3_(),
      rank_simulation_(rank_simulation)
  { reset(); }

  /// Destructor
  ~ItChild() throw()
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    TRACEPUP;
    PUParray(p,ic3_,3);
    p | rank_simulation_;
  }

  /// Go to next child, returning false when done
  bool next (int ic3[3]) throw()
  {

    increment_() ;

    ic3[0] = rank_simulation_ >= 1 ? ic3_[0] : 0;
    ic3[1] = rank_simulation_ >= 2 ? ic3_[1] : 0;
    ic3[2] = rank_simulation_ >= 3 ? ic3_[2] : 0;
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
      if (rank_simulation_ >= 1 && ic3_[0] < 1) {
	++ic3_[0];
      }	else {
	ic3_[0] = 0;
	if (rank_simulation_ >= 2 && ic3_[1] < 1) {
	  ++ic3_[1];
	} else {
	  ic3_[1] = 0;
	  if (rank_simulation_ >= 3 && ic3_[2] < 1) {
	    ++ic3_[2];
	  } else {
	    reset();
	  }
	}
      }
    }
  }

  void set_first_()
  {
    ic3_[0] = 0;
    ic3_[1] = 0;
    ic3_[2] = 0;
  }

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Current child
  int ic3_[3];

  /// simulation rank
  int rank_simulation_;

};

#endif /* MESH_ITCHILD_HPP */

