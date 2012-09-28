// See LICENSE_CELLO file for license and copyright information

/// @file     Io_IoFieldBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoFieldBlock class

#include "io.hpp"

//----------------------------------------------------------------------

IoFieldBlock::IoFieldBlock() throw ()
  : Io(4,1),
    field_descr_(0),
    field_block_(0),
    field_index_(0)

{
  meta_name_.push_back("size");
  meta_name_.push_back("array_size");
  meta_name_.push_back("num_fields");
  meta_name_.push_back("ghosts_allocated");

  ASSERT2("IoFieldBlock::IoFieldBlock()",
	 "meta_name.size() [%d] !=  meta_count_ [%d]",
	  meta_name_.size(),meta_count(),
	  meta_name_.size()==meta_count());
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
void IoFieldBlock::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change

  Io::pup(p);

  if (p.isUnpacking()) field_descr_ = new FieldDescr;
  p | *field_descr_;
  if (p.isUnpacking()) field_block_ = new FieldBlock;
  p | *field_block_;
  p | field_index_;
}
#endif

//----------------------------------------------------------------------

void IoFieldBlock::meta_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{
}

//----------------------------------------------------------------------

void IoFieldBlock::data_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
  if (buffer) (*buffer) = (void * ) 
		field_block_->field_unknowns(field_descr_,field_index_);
  if (name)   (*name)   = field_descr_->field_name(field_index_);
  if (type) {

    precision_type precision = field_descr_->precision(field_index_);
    if (precision == precision_default) precision = default_precision;

    switch (precision) {
    case precision_single: (*type) = scalar_type_float; break;
    case precision_double: (*type) = scalar_type_double; break;
    default:
      ERROR2 ("IoFieldBlock",
	      "Unsupported precision type %d for field %s",
	      precision, field_descr_->field_name(field_index_).c_str());
    }
  }
  int nbx,nby,nbz;
  field_block_->size(&nbx,&nby,&nbz);

  int ngx=0,ngy=0,ngz=0;

  field_descr_->ghosts(field_index_,&ngx,&ngy,&ngz);

  // Exclude ghosts when writing

  if (field_block_->ghosts_allocated()) {

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
