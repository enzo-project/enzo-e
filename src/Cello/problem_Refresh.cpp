// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-08-11
/// @brief    

#include "problem.hpp"

//----------------------------------------------------------------------

void Refresh::add_field(std::string field_name)
{
  const int id_field = cello::field_descr()->field_id(field_name);
  add_field(id_field);
}

//----------------------------------------------------------------------

int Refresh::data_size () const
{
  int count = 0;

  // WARNING: Skipping many fields since data methods are only called
  // when the Refresh object is a member of FieldFace, which in turn
  // only accesses field and particle lists and accumulate_

  SIZE_ARRAY(&count,field_list_src_);
  SIZE_ARRAY(&count,field_list_dst_);
  SIZE_ARRAY(&count,particle_list_);
  
  SIZE_VALUE(&count,all_fields_);
  SIZE_VALUE(&count,all_particles_);
  SIZE_VALUE(&count,accumulate_);

  return count;

}

//----------------------------------------------------------------------

char * Refresh::save_data (char * buffer) const
{
  char * p = buffer;

  SAVE_ARRAY(&p,field_list_src_);
  SAVE_ARRAY(&p,field_list_dst_);
  SAVE_ARRAY(&p,particle_list_);
  
  SAVE_VALUE(&p,all_fields_);
  SAVE_VALUE(&p,all_particles_);
  SAVE_VALUE(&p,accumulate_);

  ASSERT2 ("Refresh::save_data\n",
 	   "Actual size %d does not equal computed size %d",
	   p-buffer,data_size(),
	   ((p-buffer)==data_size()));

  return p;
  
  
}

//----------------------------------------------------------------------

char * Refresh::load_data (char * buffer)
{
  char * p = buffer;

  LOAD_ARRAY(&p,field_list_src_);
  LOAD_ARRAY(&p,field_list_dst_);
  LOAD_ARRAY(&p,particle_list_);

  LOAD_VALUE(&p,all_fields_);
  LOAD_VALUE(&p,all_particles_);
  LOAD_VALUE(&p,accumulate_);

  ASSERT2 ("Refresh::load_data\n",
	   "Actual size %d does not equal computed size %d",
	   p-buffer,data_size(),
	   ((p-buffer)==data_size()));

  return p;
}

