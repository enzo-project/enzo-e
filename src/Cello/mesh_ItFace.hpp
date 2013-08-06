// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItFace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    [\ref Mesh] Declaration of the ItFace class
///

#ifndef MESH_ITFACE_HPP
#define MESH_ITFACE_HPP

class ItFace {

  /// @class    ItFace
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  ItFace(int rank_simulation, int rank_limit) throw();

  /// Destructor
  ~ItFace() throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    TRACEPUP;
    PUParray(p,if3_,3);
    p | rank_simulation_;
    p | rank_limit_;
  }
#endif

  /// Reduce another value
  bool next (int if3[3]) throw();

  /// Reset the Iterator to the beginning
  void reset() throw();

  bool is_reset() const;

  /// Return the current value of the reduction operator
  void value(int if3[3]) const throw();

private: // functions

  /// go to the next 
  void increment_();

  void set_first_();

  /// Whether the current face rank is valid
  bool valid_() const;

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Current face
  int if3_[3];

  /// simulation rank
  int rank_simulation_;

  /// face rank limit
  int rank_limit_;

};

#endif /* MESH_ITFACE_HPP */

