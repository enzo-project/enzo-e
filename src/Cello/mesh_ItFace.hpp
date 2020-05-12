// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItFace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    [\ref Mesh] Declaration of the ItFace class
///

#ifndef MESH_IT_FACE_HPP
#define MESH_IT_FACE_HPP

class ItFace {

  /// @class    ItFace
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  ItFace(int rank,
	 int rank_limit,
	 bool periodic[3],
	 int n3[3],
	 Index index,
	 const int * ic3=0,
	 const int * if3=0) throw();

  /// Destructor
  ~ItFace() throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    TRACEPUP;
    PUParray(p,if3_,3);
    p | ic3_;
    p | ipf3_;
    p | rank_;
    p | rank_limit_;
    PUParray (p,periodicity_,3);
    PUParray (p,n3_,3);
    p | index_;
  }

  /// Go to the next face if any and return it through of3[]
  bool next (int of3[3]) throw()
  {
    const bool retval = next_();
    if (retval) face_(of3);
    return retval;
  }

  Index index() const ;

  /// Reset the Iterator to the beginning
  void reset() throw();

  bool is_reset() const;

private: // functions

  /// Go to the next face if any
  bool next_ () throw();

  /// Return the current face through of3[]
  void face_ (int of3[3]) const ;

  /// go to the next face
  void increment_();


  /// go to the first face
  void set_first_();

  /// Whether the current face rank is valid
  bool valid_() const;

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Current face
  int if3_[3];

  /// Adjacency child
  std::vector<int> ic3_;

  /// Parent face
  std::vector<int> ipf3_;

  /// simulation rank
  int rank_;

  /// face rank limit
  int rank_limit_;

  /// Periodicity
  int periodicity_[3];

  /// Size of the octree array
  int n3_[3];

  /// Index
  Index index_;

};

#endif /* MESH_ITFACE_HPP */

