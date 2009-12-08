//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Point class for representing a single point

/**
 * 
 * @file      point.cpp
 * @brief     Implementation of the Point class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 */

#include <assert.h>
#include <stdio.h>
#include <string>

#include "newgrav-scalar.h"
#include "newgrav-point.h"

//----------------------------------------------------------------------

int Point::d_ = 0;

//----------------------------------------------------------------------

Point::Point (Scalar m,Scalar x1, Scalar x2, Scalar x3) throw ()
  : m_(m),
    igrid_(0)
{ 
  assert (d_ > 0); 
  alloc_();
  if (d_ >= 1) x_[0] = x1;
  if (d_ >= 2) x_[1] = x2;
  if (d_ >= 3) x_[2] = x3;
}

//----------------------------------------------------------------------

Point::Point (std::string parms) throw ()
{
  // Define a point given text parameters, typically from a file

  alloc_();

  assert (d_ >= 3);
  input (parms);
}
//----------------------------------------------------------------------

Point::~Point () throw ()
{
  dealloc_();
}
  
//======================================================================

void Point::print () throw ()
{
  printf ("Point\n");
  printf ("   position   ");
  if (d_>=1) printf (SCALAR_PRINTF, x_[0]);
  if (d_>=2) printf (SCALAR_PRINTF, x_[1]);
  if (d_>=3) printf (SCALAR_PRINTF, x_[2]);
  printf ("\n");
  printf ("   mass       " SCALAR_PRINTF "\n",m_);
  printf ("   grid        %d\n",igrid_);
}

//----------------------------------------------------------------------

void Point::write (FILE *fp) throw ()
{
  assert (d_ == 3);
  if (fp == 0) fp = stdout;

  fprintf (fp,"point "
	   SCALAR_SCANF " "
	   SCALAR_SCANF " "SCALAR_SCANF " "SCALAR_SCANF 
	   "%d\n",
	   m_,
	   x_[0],x_[1],x_[2],
	   igrid_);
}

//----------------------------------------------------------------------

void Point::input (std::string parms) throw ()
{
  assert (d_ == 3);
  sscanf (parms.c_str(),
	  SCALAR_SCANF 
	  SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF "%d", 
	  &m_,
	  &x_[0],&x_[1],&x_[2],
	  &igrid_);
}

//======================================================================

void Point::alloc_ () throw ()
{
  assert (d_ > 0);
  x_ = new Scalar [ d_ ];
}

//----------------------------------------------------------------------

void Point::dealloc_ () throw ()
{
  delete [] x_;   x_ = 0;
}
