// $Id: amr.hpp 1259 2010-03-02 03:12:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     amr_amr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Nov 10 15:38:40 PST 2009
/// @brief    Declaration of the Amr class

#ifndef AMR_AMR_HPP
#define AMR_AMR_HPP

#include <memory>

class Amr {

  /// @class    Amr
  /// @ingroup  Amr
  /// @brief    Adaptive mesh refinement hierarchy

public: // interface

  /// Initialize an Amr object
  Amr() :
    max_level_(0),
    tree_(0)
  {};

private: // attributes

  /// Maximum level for the hierarchy (0 = unigrid)
  int max_level_;

  /// Tree defining the AMR hierarchy topology
  strict_auto_ptr<TreeK> tree_;


};

#endif /* AMR_AMR_HPP */

