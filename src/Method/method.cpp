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


/** 
 *********************************************************************
 *
 * @file      method_method.cpp
 * @brief     Implements the Method base class
 * @author    James
 * @date      Mon Jul 13 11:12:25 PDT 2009
 *
 * DESCRIPTION 
 * 
 *    Implements the Method base class
 *
 * PACKAGES
 *
 *    Field
 *    Particle
 * 
 * INCLUDES
 *  
 *    
 *
 * PUBLIC FUNCTIONS
 *  
 *    Method::Method()
 *    Method::initialize()
 *    Method::add_argument()
 *    Method::apply()
 *
 * PRIVATE FUCTIONS
 *  
 * 
 *
 * $Id$
 *
 *********************************************************************
 */

#include <string>

#include "method.hpp"
 
/**
 *********************************************************************
 *
 * @return        There is no return value
 *
 * Create a new Method
 *
 *********************************************************************
 */

Method::Method() throw()
{
}

/**
 *********************************************************************
 *
 * @param  method The Method being copied 
 * @return        There is no return value
 *
 * Create a copy of the given Method
 *
 *********************************************************************
 */

Method::Method(const Method & method) throw()
{
}

/**
 *********************************************************************
 *
 * @return        There is no return value
 *
 * Initialize the Method
 *
 *********************************************************************
 */

void Method::initialize() throw()
{
}

/**
 *********************************************************************
 *
 * @param  name   Name of the Field or Particle set, or tag
 * @return        There is no return value
 *
 * Specify a Field or Particle type read or modified 
 *
 *********************************************************************
 */

void Method::add_argument(std::string name) throw()
{
}

/**
 *********************************************************************
 *
 * @return        There is no return value
 *
 * Apply the method
 *
 *********************************************************************
 */

void Method::apply() throw()
{
}
    
