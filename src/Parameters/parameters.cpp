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
 * @file      Parameters.cpp
 * @brief     Read in a parameter file and access parameter values
 * @author    James Bordner
 * @date      Thu Jul  9 15:38:43 PDT 2009
 * @bug       
 * @note      
 *
 * DESCRIPTION 
 * 
 *    Class Parameters to read in a parameter file and access
 *    parameter values
 *
 * PACKAGES
 *
 *    NONE
 * 
 * INCLUDES
 *  
 *    parameters.hpp
 *
 * PUBLIC FUNCTIONS
 *  
 *    Parameters::Parameters()
 *    Parameters::~Parameters()
 *    Parameters::read()
 *    Parameters::set_value()
 *    Parameters::get_string()
 *
 * PRIVATE FUCTIONS
 *  
 * $Id$
 *
 *********************************************************************
 */

#include <string.h>

#include "parameters.hpp"
 
/**
 *********************************************************************
 *
 * @param   
 * @return  
 *
 * This function creates an empty Parameters object
 *
 *********************************************************************
 */

Parameters::Parameters() 
  throw()
  : values_()
{
}

/**
 *********************************************************************
 *
 * @param   file_pointer: An opened input parameter file or stdin
 * @return  There is no return value
 *
 * This function reads in parameter-value key pairs, one per line
 *
 *********************************************************************
 */

void
Parameters::read 
( FILE * file_pointer ) 
  throw(ExceptionBadPointer)

{
  // check input

  
  if (file_pointer == 0) {
    printf ("Throwing exception\n");
    throw ExceptionBadPointer();
    printf ("Threw exception??\n");
  }

  char buffer[MAX_PARAMETER_FILE_WIDTH];
  char parameter[MAX_PARAMETER_FILE_WIDTH];
  
  // Read in one line at a time
  while (readline_(file_pointer,buffer,MAX_PARAMETER_FILE_WIDTH)) {

    // Get the parameter
    sscanf(buffer,"%s",parameter);

    // Find the value
    const char * value = buffer + strlen(parameter);
    while (isspace(*++value)) ;

    add_parameter_(parameter,value);

  }
}

/**
 *********************************************************************
 *
 * @param   parameter
 * @return  Return the value of the parameter if it exists, or "" if not
 *
 * Return the value of the parameter if it exists, or "" if not
 *
 *********************************************************************
 */

std::string 
Parameters::get_string 
( std::string parameter  ) 
  throw()

{
  return (values_.find(parameter) == values_.end()) ? 
    "" : values_.find(parameter)->second;
}

/**
 *********************************************************************
 *
 * @param   file_pointer       the opened input file
 * @param   buffer             a character string buffer
 * @param   buffer_length      size of the buffer
 *
 * @return  Whether the end of the file has been reached
 *
 * Read in the next line of the input file into the buffer
 *
 *********************************************************************
 */

int 
Parameters::readline_ 
( FILE*  file_pointer, 
  char * buffer, 
  int    buffer_length ) 
  throw()

{

  int i;

  // Clear the line buffer

  for (i=0; i<buffer_length; i++) {
    buffer[i]=0;
  }

  // Read the next line into the buffer

  int c = 0;
  for (i=0; c != EOF && c != '\n' && i < buffer_length-1; i++) {
    buffer[i] = c = fgetc(file_pointer);
  }

  // Back up i to last character read
  i--;

  // Check for buffer overrun

  if (i+1 >= buffer_length-1) {
    throw "Input file line too long";
  }

  // Convert the buffer into a C string

  if (buffer[i] == '\n') buffer[i] = '\0';

  // Return whether or not the end-of-file has been reached

  return (c != EOF);

}

void    
Parameters::add_parameter_ 
( std::string parameter,
  std::string value ) 
  throw()

{
    if (parameter=="" || parameter=="#" || parameter=="//") return;
  values_.insert( std::pair<std::string,std::string>(parameter,value) );
}
