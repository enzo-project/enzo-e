#ifndef METHOD_HPP
#define METHOD_HPP

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
 * SYNOPSIS:
 *
 *    
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
 * SYNOPSIS:
 *
 *    
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
//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

#include "cello.h"

/** 
 *********************************************************************
 *
 * @file      method_method.hpp
 * @brief     Defines the Method base class
 * @author    James Bordner
 * @date      Mon Jul 13 11:11:47 PDT 2009
 *
 * Defines the Method base class
 *
 * $Id$
 *
 *********************************************************************
 */

class Method {

/** 
 *********************************************************************
 *
 * @class     Method
 * @brief     Base class for numerical methods
 * @ingroup   Method
 *
 * Base class for numerical methods
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

  /// Create a new Method
  Method() throw();

  /// Create Method
  Method(const Method &) throw();

  /// Initialize the method
  void initialize() throw(); 

  /// Specify a Field or Particle type read or modified 
  void add_argument(std::string) throw(); 

  /// Apply the method
  void apply() throw(); 

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

};

#include "method_ppm.hpp"

#endif /* METHOD_HPP */
