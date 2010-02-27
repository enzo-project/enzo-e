// $Id:
mpi.hpp 1254 2010-02-26 00:35:30Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MPI_HPP
#define MPI_HPP

/// @file
/// @brief     
/// @author    
/// @date      
///
/// Detailed description of file mpi.hpp


/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
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
 * @note      
 *
 * MPI helper functions
 *
 * $Id$
 *
 *********************************************************************
 */

class Mpi {

  /// @class    Foo
  /// @brief    Brief description of class Foo.
  /// @ingroup  Template

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

