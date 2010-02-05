#ifndef AMR_HPP
#define AMR_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      amr.hpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @package   Amr
 * @date      Tue Nov 10 15:38:40 PST 2009
 *
 * Adaptive mesh refinement hierarchy
 *
 * $Id$
 *
 *********************************************************************
 */

#include "amr_treek.hpp"

class Amr {

/** 
 *********************************************************************
 *
 * @class     Amr
 * @brief     
 * @ingroup   Amr
 *
 * 
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// 
  Amr() :
    max_level_(0),
    nx0_(0),ny0_(0),nz0_(0),
    tree_(0)
  {};

  /// 
  ~Amr() 
  {
    delete tree_;
  };

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// 

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Maximum level for the hierarchy (0 = unigrid)
  int max_level_;

  /// Size of the root grid
  int nx0_, ny0_, nz0_;

  /// Tree defining the AMR hierarchy topology
  TreeK * tree_;


};

#endif /* AMR_HPP */

