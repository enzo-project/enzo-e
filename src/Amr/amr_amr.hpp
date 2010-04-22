// $Id: amr.hpp 1259 2010-03-02 03:12:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef AMR_AMR_HPP
#define AMR_AMR_HPP

/// @file     amr_amr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Nov 10 15:38:40 PST 2009
/// @brief    Declaration of the Amr class

class Amr {

  /// @class    Amr
  /// @ingroup  Amr
  /// @brief    Adaptive mesh refinement hierarchy

public: // interface

  /// Initialize an Amr object
  Amr() :
    max_level_(0),
    nx0_(0),ny0_(0),nz0_(0),
    tree_(0)
  {};

  /// Delet an Amr object
  ~Amr() 
  { delete tree_; };

  /// Copy constructor
  Amr(const Amr & amr) throw()
  { INCOMPLETE_MESSAGE("Amr::Amr",""); };

  /// Assignment operator
  Amr & operator= (const Amr & amr) throw()
  { INCOMPLETE_MESSAGE("Amr::operator =",""); };


private: // attributes

  /// Maximum level for the hierarchy (0 = unigrid)
  int max_level_;

  /// Size of the root grid
  int nx0_, ny0_, nz0_;

  /// Tree defining the AMR hierarchy topology
  TreeK * tree_;


};

#endif /* AMR_AMR_HPP */

