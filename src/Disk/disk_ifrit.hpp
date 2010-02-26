//
// $Id$
//
// See LICENSE_CELLO file for license and copyright information
//

#ifndef DISK_IFRIT_HPP
#define DISK_IFRIT_HPP

///
/// @brief     
/// @author    
/// @date      
/// @ingroup
/// @note      
///

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
 * @file      disk_ifrit.hpp
 * @brief     Definition of the Ifrit class
 * @author    James Bordner
 * @date      Thu Feb 21 16:05:34 PST 2008
 *
 * $Id$
 *
 *********************************************************************
 */
 
class Ifrit {

/** 
 *********************************************************************
 *
 * @class     Ifrit
 * @brief     Class for writing and reading IFRIT files
 * @ingroup   Storage
 *
 * An Ifrit object currently corresponds to a single IFRIT file / group
 * / dataset.
 *
 *********************************************************************
 */

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

public:

  /// Initialize the Ifrit object
  Ifrit() {};
  ~Ifrit() {};
  void read_bin  (std::string name, Scalar * buffer, int * nx, int * ny, int * nz);
  void write_bin (std::string name, Scalar * buffer, int   nx, int   ny, int   nz);

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

private:

};

#endif /* DISK_IFRIT_HPP */

