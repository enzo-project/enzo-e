//
// $Id$
//
// See LICENSE_CELLO file for license and copyright information
//

#ifndef MPI_HPP
#define MPI_HPP

/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
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


/** 
 *********************************************************************
 *
 * @file      mpi.hpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      Thu Oct 15 10:40:37 PDT 2009
 * @bug       
 * @note      
 *
 * MPI helper functions
 *
 * $Id$
 *
 *********************************************************************
 */

class Mpi {

/** 
 *********************************************************************
 *
 * @class     Mpi
 * @brief     MPI helper functions
 * @ingroup   Parallel
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
  Mpi();

  /// 
  ~Mpi();

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// 

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// 
  bool blocking_;
  
  ///
  enum type_send {
    type_standard,
    type_buffered,
    type_synchronous,
    type_ready } type_;
  


};

#endif /* MPI_HPP */

