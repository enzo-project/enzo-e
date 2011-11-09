// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItReduce.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-26
/// @brief    Implementation of the ItReduce class
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

ItReduce * ItReduce::create (reduce_enum reduce)
{
  switch (reduce) {
  case reduce_min: return new ItReduceMin;
  case reduce_max: return new ItReduceMax;
  case reduce_avg: return new ItReduceAvg;
  case reduce_sum: return new ItReduceSum;
  default: 
    ERROR1("ItReduce::create","Unsupported reduce_enum value %d",reduce); 
    return 0;
  }
}
