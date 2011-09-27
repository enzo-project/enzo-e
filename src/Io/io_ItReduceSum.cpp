// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItReduceSum.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

ItReduceSum::ItReduceSum () throw ()
  : ItReduce(0.0)
{
}

//----------------------------------------------------------------------

ItReduceSum::~ItReduceSum () throw ()
{
}

//----------------------------------------------------------------------
 
void ItReduceSum::next (double value) throw()
{
  value_ += value;
}

//----------------------------------------------------------------------

void ItReduceSum::first() throw()
{
  value_ = 0.0;
}

//----------------------------------------------------------------------

double ItReduceSum::value() const throw()
{
  return value_;
}

//----------------------------------------------------------------------

//======================================================================

