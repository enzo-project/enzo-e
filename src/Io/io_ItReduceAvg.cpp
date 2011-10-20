// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItReduceAvg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

ItReduceAvg::ItReduceAvg () throw ()
  : ItReduce(0.0)
{
}

//----------------------------------------------------------------------

ItReduceAvg::~ItReduceAvg () throw ()
{
}

//----------------------------------------------------------------------
 
void ItReduceAvg::next (double value) throw()
{
  count_++;
  value_ += value;
}

//----------------------------------------------------------------------

void ItReduceAvg::first() throw()
{
  count_=0;
  value_ = 0.0;
}

//----------------------------------------------------------------------

double ItReduceAvg::value() const throw()
{
  return value_ / count_;
}

//----------------------------------------------------------------------

//======================================================================

