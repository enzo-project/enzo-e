// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_CountersUser.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-20
/// @brief    Implementation of the CountersUser class

#include "lcaperf.hpp"

//----------------------------------------------------------------------

namespace lca {

//----------------------------------------------------------------------

CountersUser::CountersUser() throw ()
  : Counters(0),
    value_()
{
  TRACE("CountersUser::CountersUser");
}

//----------------------------------------------------------------------

CountersUser::~CountersUser() throw ()
{
  TRACE("CountersUser::~CountersUser");
}

//----------------------------------------------------------------------

CountersUser::CountersUser(const CountersUser & counters) throw ()
  : Counters(counters)
/// @param     counters  Object being copied
{
  TRACE("CountersUser::CountersUser");
}

//----------------------------------------------------------------------

CountersUser & CountersUser::operator= (const CountersUser & counters) throw ()
/// @param     counters  Source object of the assignment
/// @return    The target assigned object
{
  TRACE("CountersUser & CountersUser::operator= ");
  return *this;
}

//======================================================================

void CountersUser::create (std::string counter, enum counter_type type)
{
  TRACE("void CountersUser::create ");
  
  int index = num_counters_;

  ++ num_counters_;

  name_.resize(num_counters_);
  value_.resize(num_counters_);
  type_.resize(num_counters_);

  index_[counter] = index;
  name_ [index]   = counter;
  value_[index]   = 0;
  type_ [index]   = type;

}

//----------------------------------------------------------------------

void CountersUser::remove (std::string counter)
{
  TRACE("void CountersUser::remove ");
  bool delete_all = (counter == "*");

  if (delete_all) {

    index_.clear();
    name_.clear();
    value_.clear();
    type_.clear();

  } else {

    if (index_.find(counter) == index_.end()) {

      fprintf (stderr, _ERROR 
	       "remove(\"%s\") called with non-existing user_counter!\n",
	       __FILE__,__LINE__,counter.c_str());
      fflush(stderr);

    } else {
      
      int index = index_[counter];

      index_.erase (index_.find(counter));
      name_.erase  (name_.begin() + index);
      value_.erase (value_.begin() + index);
      type_.erase (type_.begin() + index);

      // update indicies

      std::map<std::string,int>::iterator iter;

      for (iter = index_.begin();  iter != index_.end(); ++iter) {
	if (iter->second > index) --iter->second;
      }

    }
  }
  
}

//----------------------------------------------------------------------

void CountersUser::increment (std::string counter, long long value)
{
  TRACE("void CountersUser::increment ");
  if (index_.find(counter) == index_.end()) {
    fprintf (stderr, _ERROR "incrementing nonexistent counter %s!\n",
	    __FILE__,__LINE__,counter.c_str());
    fflush(stderr);
    return;
  }

  value_[index_[counter]] += value;
}

//----------------------------------------------------------------------

void CountersUser::assign (std::string counter, long long value)
{
  TRACE("void CountersUser::assign ");
  if (index_.find(counter) == index_.end()) {
    fprintf (stderr, _ERROR "assigning nonexistent counter %s!\n",
	    __FILE__,__LINE__,counter.c_str());
    return;
  }

  value_[index_[counter]] = value;
}

//======================================================================


long long * CountersUser::start_()
{
  TRACE("long long * CountersUser::start_");
  long long * counters = new long long [num_counters_];

  for (int k=0; k<num_counters_; k++) {
    counters[k] = value_[k];
  }

  return counters;
}

//----------------------------------------------------------------------

void CountersUser::stop_(long long * counters)
{
  TRACE("void CountersUser::stop_");
  for (int k=0; k<num_counters_; k++) {
    if (type_[k] == counter_type_relative) 
      counters[k] = value_[k] - counters[k];
    if (type_[k] == counter_type_absolute) 
      counters[k] = value_[k];
  }
}

//----------------------------------------------------------------------

void CountersUser::update_(std::string key, long long * counters)
{
  TRACE("void CountersUser::update_");

  if (global_.find(key) == global_.end() || (global_.at(key)==NULL)) {
    global_[key] = new long long [num_counters_];
    long long * g = global_[key];
    for (int i=0; i<num_counters_; i++) g[i]=0;
  }

  long long * globals = global_[key];

  for (int k=0; k<num_counters_; k++) {
    if (type_[k] == counter_type_relative) 
      globals[k] += counters[k];
    if (type_[k] == counter_type_absolute) 
      globals[k] = counters[k];
  }  

  // Delete local counters

  delete [] counters;
}

//======================================================================

}
