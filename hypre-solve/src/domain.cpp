
/// Domain class for representing a problem domain

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

#include <assert.h>

#include <string>

#include "scalar.hpp"
#include "domain.hpp"

//----------------------------------------------------------------------

Domain::Domain () throw ()
  : d_(0)
{ 
  for (int i=0; i<3; i++) xl_[i] = xu_[i] = 0.0;
}

//----------------------------------------------------------------------

Domain::Domain (std::string parms) throw ()
{
  // Define a domain given text parameters, typically from an input file

  read (parms);
  assert (1 <= d_ && d_ <= 3);
}
//----------------------------------------------------------------------

Domain::~Domain () throw ()
{
}
  
//======================================================================

void Domain::print () throw ()
{
  printf ("Domain\n"
	  "   dimension      %d\n"
	  "   lower position "SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF"\n"
	  "   upper position "SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF"\n",
	  d_,
	  xl_[0],xl_[1],xl_[2],
	  xu_[0],xu_[1],xu_[2]);
}

//----------------------------------------------------------------------

void Domain::write (FILE *fp) throw ()
{
  assert (d_ == 3);
  if (fp == 0) fp = stdout;

  fprintf (fp,"domain "
	   "%d "
	   SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF
	   SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF,
	   d_,
	   xl_[0],xl_[1],xl_[2],
	   xu_[0],xu_[1],xu_[2]);
}

//----------------------------------------------------------------------

void Domain::read (std::string parms) throw ()
{
  sscanf (parms.c_str(),
	  "%d"
	  SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF
	  SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF,
	  &d_,&xl_[0],&xl_[1],&xl_[2],&xu_[0],&xu_[1],&xu_[2]);
}
