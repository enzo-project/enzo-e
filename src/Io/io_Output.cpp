// See LICENSE_CELLO file for license and copyright information

/// @file     io_Output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      time interval scheduling gets confused if multiple outputs scheduled per timestep
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Output class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

Output::Output () throw()
  : io_simulation_(0),
    schedule_(new Schedule),
    process_write_(0),
#ifdef CONFIG_USE_CHARM
    count_reduce_(0),
#endif
    file_name_(""),
    file_vars_(),
    field_list_()
{
}

//----------------------------------------------------------------------

void Output::set_field_list (std::vector<int> field_list) throw()
{
  field_list_ = field_list;
}

//----------------------------------------------------------------------

std::string Output::expand_file_name 
(
 int cycle,
 double time
) const throw()
{
  char buffer_curr[CELLO_STRING_LENGTH];
  char buffer_next[CELLO_STRING_LENGTH];

  // copy file name (including format strings) to buffer

  strcpy (buffer_curr,file_name_.c_str());

  // loop through file_vars_[] and replace cycle or time variables
  for (size_t i=0; i<file_vars_.size(); i++) {
    if (file_vars_[i] == "cycle") {
      sprintf (buffer_next,buffer_curr, cycle);
    } else if (file_vars_[i] == "time") {
      sprintf (buffer_next,buffer_curr, time);
    } else {
      char buffer[CELLO_STRING_LENGTH];
      sprintf (buffer,"Unknown file variable #%d '%s' for file '%s'",
	       i,file_vars_[i].c_str(),file_name_.c_str());
      ERROR("Output::expand_file_name",buffer);
    }
    strcpy (buffer_curr, buffer_next);
  }
  return std::string(buffer_curr);
}

//----------------------------------------------------------------------

void Output::scheduled_write
(
 const FieldDescr * field_descr,
 Hierarchy * hierarchy, 
 int cycle, 
 double time, 
 bool root_call
 ) throw()
{
  if (schedule_->write_this_cycle(cycle, time)) {
    // Write all Hierarchy fields
    for (size_t i = 0; i<field_list_.size(); i++) {
      write (field_descr, i, hierarchy,cycle,time,root_call); 
    }
  }
}

//----------------------------------------------------------------------

void Output::scheduled_write
(
 const FieldDescr * field_descr,
 Patch * patch, 
 Hierarchy * hierarchy, 
 int cycle, 
 double time, 
 bool root_call
 ) throw()
{
  if (schedule_->write_this_cycle(cycle, time)) {
    // Write all Patch fields
    for (size_t i = 0; i<field_list_.size(); i++) {
      write (field_descr, i, patch,hierarchy,cycle,time,root_call); 
    }
  }
}

//----------------------------------------------------------------------

void Output::scheduled_write
(
 const FieldDescr * field_descr,
 Block * block, 
 Patch * patch, 
 Hierarchy * hierarchy, 
 int cycle, 
 double time, 
 bool root_call
 ) throw()
{
  if (schedule_->write_this_cycle(cycle, time)) {
    // Write all Block fields
    for (size_t i = 0; i<field_list_.size(); i++) {
      write (field_descr, i, block, patch, hierarchy,cycle,time,root_call); 
    }
  }
}

