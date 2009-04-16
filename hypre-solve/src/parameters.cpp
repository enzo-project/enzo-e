//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Parameters class for representing a set of parameters read from an input file.

/**
 * 
 * @file      parameters.cpp
 * @brief     Implementation of the Parameters class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <string>
#include <map>

#include "newgrav-hypre-solve.h"
#include "newgrav-parameters.h"

const int trace = 0;

//----------------------------------------------------------------------

/// Create a Parameters object

/** Creates an empty Parameters object. */

Parameters::Parameters () throw ()
{
  _TRACE_;
}

//----------------------------------------------------------------------

/// Delete a Parameters object

/** Deletes an empty Parameters object. */

Parameters::~Parameters () throw ()
{
  _TRACE_;
}

//----------------------------------------------------------------------

/// Read a file of key-value pairs.

/** Reads a file of key-value pairs.  Allows more than one value per key. */

void Parameters::read (std::string filename) throw ()
{

  _TRACE_;

  // Open the file and initialize the line buffer

  FILE *fp = fopen (filename.c_str(),"r");
  char buffer[BUFFER_LENGTH];
  char key[BUFFER_LENGTH];
  int i;

  // Read the file line by line

  while (readline_ (fp,buffer,BUFFER_LENGTH)) {

    // get keys in key[]
    for (i=0; i<BUFFER_LENGTH; i++) key[i]=0;
    sscanf(buffer,"%s",key);
    const char * values = buffer + strlen(key);

    // Skip over initial white space
    while (isspace(*++values))
      ;

    if (strcmp(key,"include") == 0) {

      printf ("%s:%d Including %s\n",__FILE__,__LINE__,values);
      // Recursively read included input files.  No error checking

      this->read(values);

    } else {

      add_parameter(key,values);

    }

  }

}

//----------------------------------------------------------------------

/// Print all parameters to stdout.

/** Print all parameters to stdout. */

void Parameters::print () throw ()
{
  printf ("Parameters\n");
  
  for( std::multimap<std::string, std::string>::iterator iter = values_.begin(); 
       iter != values_.end(); 
       ++iter ) {
    printf ("     Key: %s\n",iter->first.c_str());
    printf ("   Value: %s\n",iter->second.c_str());
  }
}
//----------------------------------------------------------------------

/// Associate the given value with the given key.

/** Associate the given value with the given key. */

void Parameters::add_parameter (std::string key, std::string val) throw ()
{
  _TRACE_;
  // Do nothing if key is empty or a comment ("#" or "//")
  if (key=="" || key=="#" || key=="//") return;
  values_.insert( std::pair<std::string,std::string>(key,val) );
}

//----------------------------------------------------------------------

/// Retrieve the ith value of the the given parameter.

/** Retrieve the ith value of the the given parameter. */

std::string Parameters::ith_value  (std::string key, int i) const throw ()
{
  _TRACE_;
  return "";
}

//----------------------------------------------------------------------

/// Return the multiplicity of values for the given key.

/**  Return the multiplicity of values for the given key.  May be 0. */

int Parameters::num_values  (std::string key) const throw ()
{
  _TRACE_;
  return 0;
}

//----------------------------------------------------------------------

/// Return the value for the given key.

/**  Return the first value for the given key.  Returns "" if not defined. */

std::string Parameters::value (std::string key) const throw ()
{
  return (values_.find(key) == values_.end()) ? "" : values_.find(key)->second;
}

//======================================================================

int Parameters::readline_ (FILE* fp, char * buffer, int n) throw()
{
  _TRACE_;
  int i;
  // Clear the line buffer
  for (i=0; i<BUFFER_LENGTH; i++) buffer[i]=0;
  i=0;
  int c;
  buffer[i] = c = fgetc(fp);
  while (c != EOF && c != '\n' && i < n-1) {
    ++i;
    buffer[i] = c = fgetc(fp);
  }
  if (buffer[i] == '\n') buffer[i] = '\0';
  if (i == n-1) {
    fprintf (stderr,"Line too long: %s\n",buffer);
    exit(1);
  }
  _TRACE_;
  return (c != EOF);
}

