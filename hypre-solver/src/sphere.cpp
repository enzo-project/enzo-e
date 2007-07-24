//======================================================================
//
//        File: sphere.cpp
//
//     Summary: Sphere class source file
//
// Description:
//
//      Author: James Bordner <jobordner@ucsd.edu>
//
//        Date: 2007-03-26
//
//======================================================================

#include <assert.h>
#include <stdio.h>
#include <string>

#include "scalar.hpp"
#include "sphere.hpp"

//----------------------------------------------------------------------

int Sphere::d_ = 0;

//----------------------------------------------------------------------

Sphere::Sphere () throw ()
  : m_(0),
    r_(0)
{
  alloc_();
}

//----------------------------------------------------------------------

Sphere::Sphere (std::string parms) throw ()
{
  // Define a sphere given text parameters, typically from a file

  alloc_();   // NOTE: inefficient for dimension < 3

  read (parms);
}
	  
//----------------------------------------------------------------------

Sphere::Sphere (Scalar mass,
		Scalar radius, 
		Scalar c[3],
		Scalar v[3]) throw ()
  :  m_(mass),
     r_(radius)
{
  check ();
  alloc_(d_);
  if (d_ >= 1) c_[0] = c[0];
  if (d_ >= 2) c_[1] = c[1];
  if (d_ >= 3) c_[2] = c[2];
}

//----------------------------------------------------------------------

Sphere::~Sphere () throw ()
{
  dealloc_();
  m_ = 0;
  r_ = 0;
  d_ = 0;
}

//----------------------------------------------------------------------

Sphere::Sphere (const Sphere & s) throw ()
{
  m_ = s.m_;
  r_ = s.r_;
  d_ = s.d_;
  check ();
  alloc_(d_);
  if (d_ >= 1) c_[0] = s.c_[0];
  if (d_ >= 2) c_[1] = s.c_[1];
  if (d_ >= 3) c_[2] = s.c_[2];
}

//----------------------------------------------------------------------

Sphere & Sphere::operator = (const Sphere & s) throw ()
{
  m_ = s.m_;
  r_ = s.r_;
  d_ = s.d_;
  check ();
  alloc_(d_);
  if (d_ >= 1) c_[0] = s.c_[0];
  if (d_ >= 2) c_[1] = s.c_[1];
  if (d_ >= 3) c_[2] = s.c_[2];
  return *this;
}

//======================================================================

void Sphere::print () throw ()
{
  printf ("Sphere\n" 
	  "   dimension  %d\n"
	  "   radius     "SCALAR_PRINTF "\n"
	  "   mass       "SCALAR_PRINTF "\n"
	  "   position   "SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF "\n",
	  d_,
	  r_,m_,
	  c_[0],c_[1],c_[2]);
}

//======================================================================

void Sphere::write (FILE *fp) throw ()
{
  if (fp == 0) fp = stdout;

  fprintf (fp,"sphere " 
	   SCALAR_PRINTF SCALAR_PRINTF 
	   SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF "\n",
	   m_,r_,
	   c_[0],c_[1],c_[2]);
}

//----------------------------------------------------------------------

void Sphere::read (std::string parms) throw ()
{
  sscanf (parms.c_str(),
  	  SCALAR_SCANF SCALAR_SCANF 
  	  SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF,
	  &m_,&r_,
	  &c_[0],&c_[1],&c_[2]);
}

//======================================================================

void Sphere::dealloc_ () throw ()
{
  delete [] c_;  c_ = 0;
}

//----------------------------------------------------------------------

void Sphere::alloc_ (int d) throw ()
{
  c_ = new Scalar [ d ];
}
