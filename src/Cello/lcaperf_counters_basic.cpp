// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_CountersBasic.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-20
/// @brief    Implementation of the CountersBasic class

#include "lcaperf.hpp"

//----------------------------------------------------------------------

namespace lca {

enum {
  //  index_time_begin   = 0,
  //  index_time_end     = 1,
  index_time         = 0,
  index_calls        = 1,
  num_counters_basic = 2 
};

//----------------------------------------------------------------------

CountersBasic::CountersBasic() throw ()
  : Counters(num_counters_basic),
    time_begin_(-1)
{
  name_[index_time]  = "time";
  type_[index_time]  = counter_type_relative;
  index_[name_[index_time]] = index_time;

  name_[index_calls] = "calls";
  type_[index_calls] = counter_type_relative;
  index_[name_[index_calls]]= index_calls;
}

//----------------------------------------------------------------------

CountersBasic::~CountersBasic() throw ()
{
}

//----------------------------------------------------------------------

CountersBasic::CountersBasic(const CountersBasic & counters) throw ()
  : Counters(counters)
/// @param     counters  Object being copied
{
}

//----------------------------------------------------------------------

CountersBasic & CountersBasic::operator= (const CountersBasic & counters) throw ()
/// @param     counters  Source object of the assignment
/// @return    The target assigned object
{
  return *this;
}

//======================================================================

long long * CountersBasic::start_()
{
  long long * counters = new long long [num_counters_];

  counters[index_time]  = wall_time_();
  counters[index_calls] = 0;
  DEBUG1("time = %f",1.0e-6*counters[index_time]);

  return counters;
}

//----------------------------------------------------------------------

void CountersBasic::stop_(long long * counters)
{
  counters[index_time] = wall_time_() - counters[index_time];
  DEBUG1("time = %f",1.0e-6*counters[index_time]);
}

//----------------------------------------------------------------------

void CountersBasic::update_(std::string key, long long * counters)
{

  DEBUG1("update_(%s)",key.c_str());
  // Create new global counters for key if needed

  if (global_.find(key) == global_.end()) {
    DEBUG0;
    global_[key] = new long long [num_counters_];
    long long * g = global_[key];
    for (int i=0; i<num_counters_; i++) g[i]=0;
  }

  long long * globals = global_[key];

  // Update global counters

  globals[index_time]     += counters[index_time];
  ++globals[index_calls];

  // Delete local counters

  delete [] counters;
}

//======================================================================

};

