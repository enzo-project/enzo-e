//
// $Id$
//
// See LICENSE_CELLO file for license and copyright information
//

#ifndef ERROR_EXCEPTION_HPP
#define ERROR_EXCEPTION_HPP

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
 * @file      error_exception.hpp
 * @brief     Exception class hierarchy
 * @author    James Bordner
 * @date      Sun Jul 12 13:12:09 PDT 2009
 *
 * Exception class hierarchy
 *
 * $Id$
 *
 *********************************************************************
 */

class Exception {

/** 
 *********************************************************************
 *
 * @class     Exception
 * @brief     Exception bass class
 * @ingroup   Error
 *
 * Exception bass class
 *
 *********************************************************************
 */

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------


};

class ExceptionBadPointer           : public Exception {};

class ExceptionBadArrayDeallocation : public Exception {};
class ExceptionBadArrayAllocation   : public Exception {};

class ExceptionMemoryBadDeallocate  : public Exception {};
class ExceptionMemoryBadAllocate    : public Exception {};
class ExceptionMemoryBadGroupHandle : public Exception {};

class ExceptionParametersBadType    : public Exception {};

#endif /* EXCEPTION_HPP */

