// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoFieldData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoFieldData class

#include "io.hpp"

//----------------------------------------------------------------------

IoFieldData::IoFieldData() throw ()
  : Io(4,1),
    field_descr_(0),
    field_data_(0),
    field_index_(0)

{
  meta_name_.push_back("size");
  meta_name_.push_back("array_size");
  meta_name_.push_back("num_fields");
  meta_name_.push_back("ghosts_allocated");

  ASSERT2("IoFieldData::IoFieldData()",
	 "meta_name.size() [%d] !=  meta_count_ [%d]",
	  meta_name_.size(),meta_count(),
	  meta_name_.size()==meta_count());
}

//----------------------------------------------------------------------

void IoFieldData::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change

  Io::pup(p);

  //  if (p.isUnpacking()) field_descr_ = new FieldDescr;
  //  p | *field_descr_;
  WARNING ("IoFieldData::pup","skipping field_descr_");
  
  //  if (p.isUnpacking()) field_data_ = new FieldData;
  WARNING ("IoFieldData::pup","skipping field_data_");
  //  p | *field_data_;
  p | field_index_;
}

//----------------------------------------------------------------------

void IoFieldData::meta_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{
}

//----------------------------------------------------------------------

void IoFieldData::data_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
  if (buffer) (*buffer) = (void * ) 
		field_data_->values(field_index_);
  if (name)   (*name)   = field_descr_->field_name(field_index_);
  if (type) {

    precision_type precision = field_descr_->precision(field_index_);
    if (precision == precision_default) precision = default_precision;

    switch (precision) {
    case precision_single: (*type) = scalar_type_float; break;
    case precision_double: (*type) = scalar_type_double; break;
    default:
      ERROR2 ("IoFieldData",
	      "Unsupported precision type %d for field %s",
	      precision, field_descr_->field_name(field_index_).c_str());
    }
  }
  int nbx,nby,nbz;
  field_data_->size(&nbx,&nby,&nbz);

  int ngx=0,ngy=0,ngz=0;

  field_descr_->ghosts(field_index_,&ngx,&ngy,&ngz);

  // Exclude ghosts when writing

  if (field_data_->ghosts_allocated()) {

    if (nxd) (*nxd) = nbx + 2*ngx;
    if (nyd) (*nyd) = nby + 2*ngy;
    if (nzd) (*nzd) = nbz + 2*ngz;

    if (nx) (*nx) = nbx;
    if (ny) (*ny) = nby;
    if (nz) (*nz) = nbz;

  } else {

    if (nxd) (*nxd) = nbx;
    if (nyd) (*nyd) = nby;
    if (nzd) (*nzd) = nbz;

    if (nx) (*nx) = nbx;
    if (ny) (*ny) = nby;
    if (nz) (*nz) = nbz;

  }
}
//----------------------------------------------------------------------
