// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-09-04
/// @brief    

#include "problem.hpp"

double Method::courant_global = 1.0;

//----------------------------------------------------------------------

Method::~Method() throw()
{
  delete schedule_;
  for (size_t i=0; i<refresh_list_.size(); i++) {
    delete refresh_list_[i];
    refresh_list_[i] = 0;
  }

}

//----------------------------------------------------------------------

void Method::pup (PUP::er &p)
{ TRACEPUP;
  PUP::able::pup(p);

  bool pk = p.isPacking();
  bool up = p.isUnpacking();

  int n;
  if (pk) n=refresh_list_.size();
  p | n;
  if (up) refresh_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | refresh_list_[i]; // PUP::able
  }

  p | refresh_list_;
  p | schedule_; // pupable
  p | courant_;

  p | required_fields_; // std::vector<str> required fields
}

//----------------------------------------------------------------------

void Method::set_schedule (Schedule * schedule) throw()
{ 
  if (schedule_) delete schedule_;
  schedule_ = schedule;
}

//----------------------------------------------------------------------

void Method::define_fields () throw()
{
  /* Ensure required fields are defined for this method */

  FieldDescr * field_descr = cello::field_descr();
  Config   * config  = (Config *) cello::config();

  bool added_fields = false; 

  for (int ifield = 0; ifield < required_fields_.size(); ifield++){
    if( ! field_descr->is_field( required_fields_[ifield] )){
      field_descr->insert_permanent( required_fields_[ifield] );

      field_descr->set_precision(ifield, config->field_precision);

      added_fields = true;
    }
  }

  // Need to reconstruct history if new fields added
  if (added_fields) field_descr->reset_history(config->field_history);
}

//----------------------------------------------------------------------

void Method::define_group_fields (std::vector<std::string> group_fields,
                                  std::string groupname) throw()
{
  /* Ensure fields are grouped correctly */

  FieldDescr * field_descr = cello::field_descr();

  for (int ifield = 0; ifield < group_fields.size(); ifield++){

    // Maybe just throw error here to keep this fully separate from above
    if( ! field_descr->is_field( required_fields_[ifield] )){
      field_descr->insert_permanent( required_fields_[ifield] );
    }

    if (!(field_descr->groups()->is_in( group_fields[ifield], groupname)) ){
      field_descr->groups()->add( group_fields[ifield], groupname);
    }

  }


}

//======================================================================

