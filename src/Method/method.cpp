/// @file      method_method.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Jul 13 11:12:25 PDT 2009
/// @brief     Implements the Method base class

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
    
