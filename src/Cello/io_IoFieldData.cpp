// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoFieldData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoFieldData class

#include "io.hpp"

//----------------------------------------------------------------------

IoFieldData::IoFieldData() throw ()
  : Io(1),
    field_data_(0),
    field_index_(0)

{
  meta_name_.push_back("size");
  meta_name_.push_back("array_size");
  meta_name_.push_back("num_fields");
  meta_name_.push_back("ghosts_allocated");
}

//----------------------------------------------------------------------

void IoFieldData::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change

  Io::pup(p);

  //  if (p.isUnpacking()) field_data_ = new FieldData;
  WARNING ("IoFieldData::pup","skipping field_data_");
  //  p | *field_data_;
  p | field_index_;
}

//----------------------------------------------------------------------

void IoFieldData::meta_value
(int index,
 void ** buffer, std::string * name, int * type,
 int * nxd, int * nyd, int * nzd) throw()
{
}

//----------------------------------------------------------------------

void IoFieldData::field_array
(int index, // WARNING: index ignored
 void ** buffer, std::string * name, int * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
  FieldDescr * field_descr = cello::field_descr();
  
  if (buffer) (*buffer) = (void * ) 
		field_data_->values(field_descr,field_index_);
  if (name)   (*name) = 
		std::string("field_") +	field_descr->field_name(field_index_);
  if (type) {

    precision_type precision = field_descr->precision(field_index_);
    if (precision == precision_default) precision = default_precision;

    switch (precision) {
    case precision_single: (*type) = type_float; break;
    case precision_double: (*type) = type_double; break;
    case precision_quadruple: (*type) = type_quadruple; break;
    default:
      ERROR2 ("IoFieldData",
	      "Unsupported precision type %d for field %s",
	      precision, field_descr->field_name(field_index_).c_str());
    }
  }
  int nbx,nby,nbz;
  field_data_->size(&nbx,&nby,&nbz);

  int ngx=0,ngy=0,ngz=0;

  field_descr->ghost_depth(field_index_,&ngx,&ngy,&ngz);

  int cx=0,cy=0,cz=0;
  field_descr->centering(field_index_,&cx,&cy,&cz);

  if (field_data_->ghosts_allocated()) {

    if (nxd) (*nxd) = nbx + 2*ngx + cx;
    if (nyd) (*nyd) = nby + 2*ngy + cy;
    if (nzd) (*nzd) = nbz + 2*ngz + cz;

    // Exclude ghosts when writing

    //    if (nx) (*nx) = nbx;
    //    if (ny) (*ny) = nby;
    //    if (nz) (*nz) = nbz;

    // Include ghosts when writing

     if (nx) (*nx) = nbx + 2*ngx + cx;
     if (ny) (*ny) = nby + 2*ngy + cy;
     if (nz) (*nz) = nbz + 2*ngz + cz;

  } else {

    if (nxd) (*nxd) = nbx + cx;
    if (nyd) (*nyd) = nby + cy;
    if (nzd) (*nzd) = nbz + cz;

    if (nx) (*nx) = nbx + cx;
    if (ny) (*ny) = nby + cy;
    if (nz) (*nz) = nbz + cz;

  }
}
//----------------------------------------------------------------------

void IoFieldData::particle_array 
(int it, int ib, int ia,
 void ** buffer, std::string * name, int * type,
 int * n, int * k) throw()
{
}

//----------------------------------------------------------------------

