#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

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
 * @file      parameters.hpp
 * @brief     Read in a parameter file and access parameter values
 * @author    James Bordner
 * @date      Thu Jul  9 15:44:21 PDT 2009
 * @bug       
 * @note      
 *
 * Class Parameters to read in a parameter file and access
 * parameter values
 *
 * $Id$
 *
 *********************************************************************
 */

#include <map>
#include <string>

#include <stdio.h>

// Maximum allowed width of a line in a parameter file
#define MAX_PARAMETER_FILE_WIDTH 255

class Parameters {

/** 
 *********************************************************************
 *
 * @class     Parameters
 * @brief     Read in a parameter file and access parameter values
 * @ingroup   Parameters
 *
 * Read in a parameter file and access parameter values
 *
 *********************************************************************
 */

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// A private attribute
  std::map <std::string,std::string> values_;

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// This function creates an empty parameters object
  Parameters() throw();

  /// This function reads in parameters from a file
  void read (FILE * file_pointer) throw();

  /// This function reads in parameters from a file
  std::string get_value (std::string parameter) throw();

  /// This function reads in parameters from a file
  void set_value (std::string parameter, std::string value) throw();

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// Read in the next line of the input file
  int readline_ (FILE* fp, char * buffer, int n) throw();

  /// Add a parameter / value pair
  void add_parameter_ ( std::string parameter,  std::string value )   throw();

};

#endif /* PARAMETERS_HPP */

