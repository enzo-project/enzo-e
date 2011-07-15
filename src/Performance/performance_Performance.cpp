// $Id: performance_Performance.cpp 2130 2011-03-20 01:00:25Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file      performance_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2009-10-15
/// @brief     Implementation of the Performance class

#include "cello.hpp"

#include "performance.hpp"

Performance::Performance ()
  : counters_(),
    attribute_names_           (NULL),
    counter_names_             (NULL),
    group_names_               (NULL),
    region_names_              (NULL),
    attribute_monotonic_       (NULL),
    current_group_             (0),
    current_region_            (0)
{
  // Create initial Counters object

  counters_.push_back(new Counters(num_attributes,num_counters));
}

//----------------------------------------------------------------------

Performance::~Performance()
{
  deallocate_();
}

//----------------------------------------------------------------------

Performance::Performance(const Performance & classname) throw()
{
  INCOMPLETE("Performance::Performance");
}

//----------------------------------------------------------------------

Performance & Performance::operator= (const Performance & classname) throw()
{
  INCOMPLETE("Performance::operator=");
  return *this;
}

//----------------------------------------------------------------------

void Performance::start () throw ()
{
  timer.start();
  papi.start();
}

//----------------------------------------------------------------------

void Performance::stop () throw ()
{
  timer.stop();
  papi.stop();
}

//----------------------------------------------------------------------

void Performance::print (const Monitor * monitor) const throw ()
{
  timer.print();
  papi.print();
  print_rusage_(monitor);
}

//----------------------------------------------------------------------

unsigned Performance::new_attribute(std::string attribute_name,
				    bool is_monotonic)
/// @param    id_attribute
/// @param    attribute_name
/// @param    type
{
  attribute_names_.push_back(attribute_name);
  attribute_monotonic_.push_back(is_monotonic);

  return attribute_names_.size()-1;
}

//----------------------------------------------------------------------

int Performance::attribute(unsigned id_attribute)
{
  INCOMPLETE("Performance::attribute");
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_attribute(unsigned id_attribute,
				int value)
{
  INCOMPLETE("Performance::set_attribute");
}

//----------------------------------------------------------------------

unsigned Performance::new_group(std::string group_name)
{
  group_names_.push_back(group_name);
  return group_names_.size()-1;
}

//----------------------------------------------------------------------

int Performance::group(unsigned id_group)
{
  INCOMPLETE("Performance::group");
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_group(unsigned id_group)
{
  INCOMPLETE("Performance::set_group");
}

//----------------------------------------------------------------------

void Performance::begin_group(unsigned group_id)
{
  
  if ( current_group_ ) {
    // begin_group() called when another group is already active
    char message [ ERROR_LENGTH ];
    sprintf (message, 
	     "Mismatch between begin_group(%s) and begin_group(%s)",
	     group_names_.at(current_group_).c_str(),
	     group_names_.at(group_id).c_str());
    WARNING("Performance::begin_group",message);

    // End the mistakenly active group
    end_group(current_group_);

  }

  current_group_ = group_id;

}

//----------------------------------------------------------------------

void Performance::end_group(unsigned id_group)
{
  if (id_group != current_group_) {
    // end_group() called with an inactive one
    char message [ ERROR_LENGTH ];
    sprintf (message, "Mismatch between begin_group(%s) and end_group(%s)",
	     group_names_[current_group_].c_str(),
	     group_names_[id_group].c_str());
    WARNING("Performance::end_group",message);
  }

  current_group_ = 0;

}

//----------------------------------------------------------------------

unsigned Performance::new_region(std::string region_name)
{
  region_names_.push_back(region_name);
  return region_names_.size()-1;
}

//----------------------------------------------------------------------

int Performance::region(unsigned id_region)
{
  INCOMPLETE("Performance::region");
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_region(unsigned id_region)
{
  INCOMPLETE("Performance::set_region");
}

//----------------------------------------------------------------------

void Performance::start_region(unsigned region_id)
{
  INCOMPLETE("Performance::start_region");
}

//----------------------------------------------------------------------

void Performance::stop_region(unsigned region_id)
{
  INCOMPLETE("Performance::stop_region");
}

//----------------------------------------------------------------------

unsigned Performance::new_counter(std::string counter_name)
{
  counter_names_.push_back(counter_name);
  return counter_names_.size()-1;
}

//----------------------------------------------------------------------

type_counter Performance::counter(unsigned id_counter)
{
  INCOMPLETE("Performance::counter");
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_counter(unsigned          id_counter,
			      type_counter value)
{
  INCOMPLETE("Performance::set_counter");
}

//----------------------------------------------------------------------

void Performance::increment_counter(unsigned          id_counter,
				    type_counter value)
{
  INCOMPLETE("Performance::increment_counter");
}

//----------------------------------------------------------------------

void Performance::flush()
{
  INCOMPLETE("Performance::flush");
}


//======================================================================
void Performance::deallocate_() throw()
{
  for (unsigned i=0; i<counters_.size(); i++) {
    delete counters_.at(i);
  }
}

//----------------------------------------------------------------------

void Performance::print_rusage_(const Monitor * monitor) const throw()
{
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);

  monitor->print ("[Performance] utime = %f",
   	  r.ru_utime.tv_sec + 
	  r.ru_utime.tv_usec * 1e-6);
  monitor->print ("[Performance] stime = %f",
   	  r.ru_stime.tv_sec + 
	  r.ru_stime.tv_usec * 1e-6);

  if (r.ru_maxrss) monitor->print ("[Performance]  maximum resident set size: %ld",  
			   1024 * r.ru_maxrss);
  if (r.ru_ixrss) monitor->print ("[Performance]  integral shared memory size: %ld",  r.ru_ixrss);
  if (r.ru_idrss) monitor->print ("[Performance]  integral unshared data size: %ld",  r.ru_idrss);
  if (r.ru_isrss) monitor->print ("[Performance]  integral unshared stack size: %ld",  r.ru_isrss);
  if (r.ru_minflt) monitor->print ("[Performance]  page reclaims (soft page faults): %ld",  r.ru_minflt);
  if (r.ru_majflt) monitor->print ("[Performance]  page faults (hard page faults): %ld",  r.ru_majflt);
  if (r.ru_nswap) monitor->print ("[Performance]  swaps: %ld",  r.ru_nswap);
  if (r.ru_inblock) monitor->print ("[Performance]  block input operations: %ld",  r.ru_inblock);
  if (r.ru_oublock) monitor->print ("[Performance]  block output operations: %ld",  r.ru_oublock);
  if (r.ru_msgsnd) monitor->print ("[Performance]  IPC messages sent: %ld",  r.ru_msgsnd);
  if (r.ru_msgrcv) monitor->print ("[Performance]  IPC messages received: %ld",  r.ru_msgrcv);
  if (r.ru_nsignals) monitor->print ("[Performance]  signals received: %ld",  r.ru_nsignals);
  if (r.ru_nvcsw) monitor->print ("[Performance]  voluntary context switches: %ld",  r.ru_nvcsw);
  if (r.ru_nivcsw) monitor->print ("[Performance]  involuntary context switches: %ld",  r.ru_nivcsw);


  monitor->print ("[Performance] hostid = %ld",gethostid());

}
