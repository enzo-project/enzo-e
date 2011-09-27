// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItReduceMin.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

ItReduceMin::ItReduceMin () throw ()
  : ItReduce(std::numeric_limits<double>::max())
{
}

//----------------------------------------------------------------------

ItReduceMin::~ItReduceMin () throw ()
{
}

//----------------------------------------------------------------------
 
void ItReduceMin::next (double value) throw()
{
  value_ = MIN(value,value_);
}

//----------------------------------------------------------------------

void ItReduceMin::first() throw()
{
  value_ = std::numeric_limits<double>::max();
}

//----------------------------------------------------------------------

double ItReduceMin::value() const throw()
{
  return value_;
}

//----------------------------------------------------------------------

//======================================================================

