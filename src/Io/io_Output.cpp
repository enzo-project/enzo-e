// See LICENSE_CELLO file for license and copyright information

/// @file     io_Output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      time interval scheduling gets confused if multiple outputs scheduled per timestep
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Output class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

Output::Output (Simulation * simulation) throw()
  : simulation_(simulation),
    file_(0),           // Initialization deferred
    schedule_(new Schedule),
    process_stride_(1), // default one file per process
#ifdef CONFIG_USE_CHARM
    count_reduce_(0),
#endif
    process_(0),        // initialization below
    cycle_(0),
    count_(0),
    time_(0),
    file_name_(""),     // set_filename()
    file_args_(),       // set_filename()
    it_field_(0)        // set_it_field()
{

  GroupProcess * group_process = GroupProcess::create();
  process_  = group_process->rank();
  delete group_process;
}

//----------------------------------------------------------------------

void Output::set_filename (std::string filename,
			   std::vector<std::string> fileargs) throw()
{
  file_name_ = filename; 
  file_args_ = fileargs;
}

//----------------------------------------------------------------------

std::string Output::expand_file_name () const throw()
{
  char buffer_curr[CELLO_STRING_LENGTH];
  char buffer_next[CELLO_STRING_LENGTH];

  // copy file name (including format strings) to buffer

  strcpy (buffer_curr,file_name_.c_str());

  // loop through file_args_[] and replace cycle or time variables
  for (size_t i=0; i<file_args_.size(); i++) {
    if (file_args_[i] == "cycle") {
      sprintf (buffer_next,buffer_curr, cycle_);
    } else if (file_args_[i] == "time") {
      sprintf (buffer_next,buffer_curr, time_);
    } else if (file_args_[i] == "count") {
      sprintf (buffer_next,buffer_curr, count_);
    } else if (file_args_[i] == "proc") {
      sprintf (buffer_next,buffer_curr, process_);
    } else {
      char buffer[CELLO_STRING_LENGTH];
      sprintf (buffer,"Unknown file variable #%d '%s' for file '%s'",
	       int(i),file_args_[i].c_str(),file_name_.c_str());
      ERROR("Output::expand_file_name",buffer);
    }
    strcpy (buffer_curr, buffer_next);
  }
  return std::string(buffer_curr);
}

//----------------------------------------------------------------------

bool Output::is_scheduled (int cycle, double time)
{
  cycle_ = cycle;
  time_  = time;
  return schedule()->write_this_cycle(cycle_,time_);
}

//----------------------------------------------------------------------

double Output::update_timestep (double time, double dt) const
{
  return schedule_->update_timestep(time,dt); 
}

//----------------------------------------------------------------------
