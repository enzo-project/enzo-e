// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-08-11
/// @brief    

#include "problem.hpp"

// #define DEBUG_REFRESH

//----------------------------------------------------------------------

void Refresh::add_field(std::string field_name)
{
  const int id_field = cello::field_descr()->field_id(field_name);
  if (id_field >= 0) {
    field_list_src_.push_back(id_field);
    field_list_dst_.push_back(id_field);
  }
}

//----------------------------------------------------------------------

int Refresh::data_size () const
{
  int count = 0;

  // WARNING: Skipping many fields since data methods are only called
  // when the Refresh object is a member of FieldFace, which in turn
  // only accesses field and particle lists and accumulate_

  count += sizeof(int);  // field_list_src_.size()
  count += sizeof(int)*field_list_src_.size();
  
  count += sizeof(int);  // field_list_dst_.size()
  count += sizeof(int)*field_list_dst_.size();
  
  count += sizeof(int);  // particle_list_.size()
  count += sizeof(int)*particle_list_.size();
  
  count += sizeof(bool); // all_fields_
  count += sizeof(bool); // all_particles_
  count += sizeof(bool); // accumulate_

  return count;

}

//----------------------------------------------------------------------

char * Refresh::save_data (char * buffer) const
{
#ifdef DEBUG_REFRESH  
  CkPrintf ("%d DEBUG_REFRESH Refresh::save_data() %p\n",
	    CkMyPe(),this);
  fflush(stdout);
#endif  
  char * p = buffer;
  int n;

  int length;

#ifdef DEBUG_REFRESH  
  CkPrintf ("DEBUG_REFRESH save_data lengths %d %d %d\n",
	    field_list_src_.size(),
	    field_list_dst_.size(),
	    particle_list_.size());
  fflush(stdout);
#endif  

  length = field_list_src_.size();
  memcpy(p,&length, n=sizeof(int));
  p+=n;
  memcpy(p,&field_list_src_[0],n=length*sizeof(int));
  p+=n;

  length = field_list_dst_.size();
  memcpy(p,&length, n=sizeof(int));
  p+=n;
  memcpy(p,&field_list_dst_[0],n=length*sizeof(int));
  p+=n;
  
  length = particle_list_.size();
  memcpy(p,&length, n=sizeof(int));
  p+=n;
  memcpy(p,&particle_list_[0],n=length*sizeof(int));
  p+=n;

#ifdef DEBUG_REFRESH
  if (field_list_src_.size()>0) {
    CkPrintf ("DEBUG_REFRESH save_data field_list_src[0] %d\n",
	      field_list_src_[0]);
  fflush(stdout);
  }
#endif  
  
  memcpy(p,&all_fields_, n=sizeof(bool));
  p+=n;
  memcpy(p,&all_particles_,n=sizeof(bool));
  p+=n;
  memcpy(p,&accumulate_,n=sizeof(bool));
  p+=n;

  ASSERT2 ("Refresh::save_data\n",
	   "Actual size %d does not equal computed size %d",
	   p-buffer,data_size(),
	   ((p-buffer)==data_size()));
  
  return p;
}

//----------------------------------------------------------------------

char * Refresh::load_data (char * buffer)
{
#ifdef DEBUG_REFRESH  
  CkPrintf ("%d DEBUG_REFRESH Refresh::load_data() %p\n",
	    CkMyPe(),this);
  fflush(stdout);
#endif  
  char * p = buffer;
  int n;

  int length;

#ifdef DEBUG_REFRESH  
  CkPrintf ("DEBUG_REFRESH load_data lengths %d %d %d\n",
	    field_list_src_.size(),
	    field_list_dst_.size(),
	    particle_list_.size());
  fflush(stdout);
#endif  

  memcpy(&length,p, n=sizeof(int));
  p+=n;
  field_list_src_.resize(length);
  memcpy(&field_list_src_[0],p,n=length*sizeof(int));
  p+=n;

  memcpy(&length,p, n=sizeof(int));
  p+=n;
  field_list_dst_.resize(length);
  memcpy(&field_list_dst_[0],p,n=length*sizeof(int));
  p+=n;
  
  memcpy(&length, p,n=sizeof(int));
  p+=n;
  particle_list_.resize(length);
  memcpy(&particle_list_[0],p,n=length*sizeof(int));
  p+=n;
  

#ifdef DEBUG_REFRESH
  if (field_list_src_.size()>0) {
    CkPrintf ("DEBUG_REFRESH save_data field_list_src[0] %d\n",
	      field_list_src_[0]);
  fflush(stdout);
  }
#endif  

  memcpy(&all_fields_,p,n=sizeof(bool));
  p+=n;
  memcpy(&all_particles_,p,n=sizeof(bool));
  p+=n;
  memcpy(&accumulate_,p,n=sizeof(bool));
  p+=n;
  
  ASSERT2 ("Refresh::load_data\n",
	   "Actual size %d does not equal computed size %d",
	   p-buffer,data_size(),
	   ((p-buffer)==data_size()));

  return p;
}

