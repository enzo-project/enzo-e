//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2009 James Bordner
 * Copyright (C) 2009 Laboratory for Computational Astrophysics
 * Copyright (C) 2009 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */


/** 
 *********************************************************************
 *
 * @file      memory.cpp
 * @brief     Functions for dynamic memory management
 * @author    James Bordner
 * @date      Thu Sep  3 16:44:18 PDT 2009
 * @bug       
 * @note      
 *
 * DESCRIPTION 
 * 
 *    Functions for dynamic memory management
 *
 * PACKAGES
 *
 *    NONE
 * 
 * INCLUDES
 *  
 *    NONE
 *
 * PUBLIC FUNCTIONS
 *  
 *   ( ) static Memory();
 *   ( ) static void * allocate(size_t size);
 *   ( ) static void * allocate(size_t size, std::string class);
 *   ( ) static deallocate();
 *   ( ) static deallocate(std::string class);
 *   ( ) static current(class);
 *   ( ) static available();
 *   ( ) static efficiency();
 *   ( ) static highest(class);
 *   ( ) static set_highest();
 *
 *
 * PRIVATE FUCTIONS
 *  
 *    
 *
 * $Id$
 *
 *********************************************************************
 */

#include "error.hpp"

#include "memory.hpp"

Memory memory;
 
Memory::Memory()
/**
 *********************************************************************
 *
 * @param         There are no parameters
 * @return        There is no return value
 *
 * Initialize the Memory class
 *
 *********************************************************************
 */
{
}

void * Memory::allocate(size_t size)
/**
 *********************************************************************
 *
 * @param  size   Number of bytes to allocate
 * @return        Pointer to the allocated memory
 *
 * Allocate memory with the default group
 *
 *********************************************************************
 */
{
  char * buffer = new char [size];
  
  return (void *) buffer;
}

void * Memory::allocate(unsigned size, std::string group)
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Memory::allocate(size,group)");
  char * buffer = new char [size];
  
  return (void *) buffer;
}

void Memory::deallocate()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Memory::deallocate()");
}

void Memory::deallocate(std::string group)
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Memory::deallocate(group)");
}

void Memory::current(std::string group)
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Memory::current(group)");
}

void Memory::available()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Memory::available()");
}

void Memory::efficiency()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Memory::efficiency()");
}

void Memory::highest(std::string group)
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Memory::highest(group)");
}

void Memory::set_highest()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  INCOMPLETE_MESSAGE("Memory::set_highest()");
}
    
