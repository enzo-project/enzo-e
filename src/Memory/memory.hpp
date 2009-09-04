#ifndef MEMORY_HPP
#define MEMORY_HPP

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
 * @file      memory.hpp
 * @brief     Functions for dynamic memory management
 * @author    James Bordner
 * @date      Thu Sep  3 16:29:56 PDT 2009
 * @bug       
 * @note      
 *
 * Functions for dynamic memory management
 *
 * $Id$
 *
 *********************************************************************
 */

#include <string>

class Memory {

/** 
 *********************************************************************
 *
 * @class     Memory
 * @brief     Brief description of the class
 * @ingroup   Group 
 *
 * Detailed description of the class
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  // Initialize the memory component
  Memory();

  //  	 Allocate memory
  void * allocate(unsigned size);
  void * allocate(unsigned size, std::string group);

  //  	De-allocate memory
  void deallocate();
  void deallocate(std::string group);

  //  	Current number of bytes allocated
  void current(std::string group);
  
  //  	Estimate of amount of local memory availables);
  void available();

  //  	Estimate of used / available memory
  void efficiency();

  //  	Maximum number of bytes allocated
  void highest(std::string group);

  //  	Specify the maximum number of bytes to use;
  void set_highest();

 private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------


private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// A private attribute
  int private_attribute_;


};

extern Memory memory;

#endif /* MEMORY_HPP */

