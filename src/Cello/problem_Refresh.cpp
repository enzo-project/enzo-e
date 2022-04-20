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

void Refresh::add_field_src_dst(std::string field_src, std::string field_dst)
{
  const int id_field_src = cello::field_descr()->field_id(field_src);
  const int id_field_dst = cello::field_descr()->field_id(field_dst);
  add_field_src_dst(id_field_src,id_field_dst);
}

//----------------------------------------------------------------------

void Refresh::add_all_fields(std::string field_group)
{
  if (field_group == "") {
    all_fields_ = true;
  } else {
    Grouping * groups = cello::field_groups();
    int n = groups->size(field_group);
    for (int i=0; i<n; i++) {
      std::string field = groups->item(field_group,i);
      add_field(field);
    }
  }
}

//----------------------------------------------------------------------

std::vector<int> Refresh::field_list_src() const
{
  std::vector<int> field_list;
  if (all_fields_) {
    int nf = cello::field_descr()->field_count();
    for (int i=0; i<nf; i++) {
      field_list.push_back(i);
    }
    return field_list;
  } else {
    return field_list_src_;
  }
}

//----------------------------------------------------------------------

std::vector<int> Refresh::field_list_dst() const
{
  std::vector<int> field_list;
  if (all_fields_) {
    int nf = cello::field_descr()->field_count();
    for (int i=0; i<nf; i++) {
      field_list.push_back(i);
    }
    return field_list;
  } else {
    return field_list_dst_;
  }
}

//----------------------------------------------------------------------

void Refresh::box_accumulate_adjust
(Box * box, int if3[3], int g3[3])
{
  if (accumulate_) {
    int gs3[3] = {0};
    for (int i=0; i<cello::rank(); i++) {
      gs3[i] = (if3[i] != 0) ? g3[i] : 0;
    }
    box->set_send_ghosts(gs3);
  }
}

//----------------------------------------------------------------------

int Refresh::coarse_padding(const Prolong * prolong) const
{
  const int pad = accumulate_ ? 0 : prolong->coarse_padding_();
  return pad;
}

//----------------------------------------------------------------------

Prolong * Refresh::prolong ()
{
  Problem * problem = cello::problem();
  Prolong * prolong = problem ? problem->prolong(id_prolong_) : nullptr;
  if (prolong == nullptr) {
    WARNING("Refresh::prolong()",
            "Creating new ProlongLinear, potential memory leak");
  }
  return prolong ? prolong : new ProlongLinear;
}

//----------------------------------------------------------------------

Restrict * Refresh::restrict ()
{
  Problem * problem = cello::problem();
  Restrict * restrict = problem ? problem->restrict(id_restrict_) : nullptr;
  if (restrict == nullptr) {
    WARNING("Refresh::restrict()",
            "Creating new RestrictLinear, potential memory leak");
  }
  return restrict ? restrict : new RestrictLinear;
}
  
//----------------------------------------------------------------------

int Refresh::data_size () const
{
  int count = 0;

  // WARNING: Skipping many fields since data methods are only called
  // when the Refresh object is a member of FieldFace, which in turn
  // only accesses field and particle lists and accumulate_

  SIZE_SCALAR_TYPE(count,int,all_fields_);
  SIZE_VECTOR_TYPE(count,int,field_list_src_);
  SIZE_VECTOR_TYPE(count,int,field_list_dst_);

  SIZE_SCALAR_TYPE(count,int,all_particles_);
  SIZE_VECTOR_TYPE(count,int,particle_list_);

  SIZE_SCALAR_TYPE(count,int,all_fluxes_);
  SIZE_SCALAR_TYPE(count,int,ghost_depth_);
  SIZE_SCALAR_TYPE(count,int,min_face_rank_);
  SIZE_SCALAR_TYPE(count,int,neighbor_type_);
  SIZE_SCALAR_TYPE(count,int,accumulate_);
  SIZE_SCALAR_TYPE(count,int,sync_type_);
  SIZE_SCALAR_TYPE(count,int,sync_id_);
  SIZE_SCALAR_TYPE(count,int,active_);
  SIZE_SCALAR_TYPE(count,int,callback_);
  SIZE_SCALAR_TYPE(count,int,root_level_);
  SIZE_SCALAR_TYPE(count,int,id_refresh_);

  SIZE_SCALAR_TYPE(count,int,id_prolong_);
  SIZE_SCALAR_TYPE(count,int,id_restrict_);

  return count;

}

//----------------------------------------------------------------------

char * Refresh::save_data (char * buffer) const
{
  char * p = buffer;

  SAVE_SCALAR_TYPE(p,int,all_fields_);
  SAVE_VECTOR_TYPE(p,int,field_list_src_);
  SAVE_VECTOR_TYPE(p,int,field_list_dst_);

  SAVE_SCALAR_TYPE(p,int,all_particles_);
  SAVE_VECTOR_TYPE(p,int,particle_list_);

  SAVE_SCALAR_TYPE(p,int,all_fluxes_);
  SAVE_SCALAR_TYPE(p,int,ghost_depth_);
  SAVE_SCALAR_TYPE(p,int,min_face_rank_);
  SAVE_SCALAR_TYPE(p,int,neighbor_type_);
  SAVE_SCALAR_TYPE(p,int,accumulate_);
  SAVE_SCALAR_TYPE(p,int,sync_type_);
  SAVE_SCALAR_TYPE(p,int,sync_id_);
  SAVE_SCALAR_TYPE(p,int,active_);
  SAVE_SCALAR_TYPE(p,int,callback_);
  SAVE_SCALAR_TYPE(p,int,root_level_);
  SAVE_SCALAR_TYPE(p,int,id_refresh_);

  SAVE_SCALAR_TYPE(p,int,id_prolong_);
  SAVE_SCALAR_TYPE(p,int,id_restrict_);

  ASSERT2 ("Refresh::save_data\n",
 	   "Actual size %ld does not equal computed size %d",
	   p-buffer,data_size(),
	   ((p-buffer)==data_size()));

  return p;
}

//----------------------------------------------------------------------

char * Refresh::load_data (char * buffer)
{
  char * p = buffer;

  LOAD_SCALAR_TYPE(p,int,all_fields_);
  LOAD_VECTOR_TYPE(p,int,field_list_src_);
  LOAD_VECTOR_TYPE(p,int,field_list_dst_);

  LOAD_SCALAR_TYPE(p,int,all_particles_);
  LOAD_VECTOR_TYPE(p,int,particle_list_);

  LOAD_SCALAR_TYPE(p,int,all_fluxes_);
  LOAD_SCALAR_TYPE(p,int,ghost_depth_);
  LOAD_SCALAR_TYPE(p,int,min_face_rank_);
  LOAD_SCALAR_TYPE(p,int,neighbor_type_);
  LOAD_SCALAR_TYPE(p,int,accumulate_);
  LOAD_SCALAR_TYPE(p,int,sync_type_);
  LOAD_SCALAR_TYPE(p,int,sync_id_);
  LOAD_SCALAR_TYPE(p,int,active_);
  LOAD_SCALAR_TYPE(p,int,callback_);
  LOAD_SCALAR_TYPE(p,int,root_level_);
  LOAD_SCALAR_TYPE(p,int,id_refresh_);

  LOAD_SCALAR_TYPE(p,int,id_prolong_);
  LOAD_SCALAR_TYPE(p,int,id_restrict_);

  ASSERT2 ("Refresh::load_data\n",
	   "Actual size %ld does not equal computed size %d",
	   p-buffer,data_size(),
	   ((p-buffer)==data_size()));

  return p;
}

