/** 
 *********************************************************************
 *
 * @file      amr.hpp
 * @brief     Adaptive mesh refinement hierarchy
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      Tue Nov 10 15:38:40 PST 2009
 * @ingroup   Amr
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */

#ifndef AMR_HPP
#define AMR_HPP

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

