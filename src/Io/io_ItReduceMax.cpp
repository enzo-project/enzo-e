// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItReduceMax.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

ItReduceMax::ItReduceMax () throw ()
  : ItReduce(std::numeric_limits<double>::min())
{
}

//----------------------------------------------------------------------

ItReduceMax::~ItReduceMax () throw ()
{
}

//----------------------------------------------------------------------
 
void ItReduceMax::next (double value) throw()
{
  value_ = MAX(value,value_);
}

//----------------------------------------------------------------------

void ItReduceMax::first() throw()
{
  value_ = std::numeric_limits<double>::min();
}

//----------------------------------------------------------------------

double ItReduceMax::value() const throw()
{
  return value_;
}

//----------------------------------------------------------------------

//======================================================================

