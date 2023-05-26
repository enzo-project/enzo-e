// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoFieldData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoFieldData class

#include "io.hpp"

//----------------------------------------------------------------------

IoFieldData::IoFieldData() throw ()
  : Io(),
    field_data_(0),
    field_index_(0),
    include_ghosts_(true)
{
  // meta_name_.push_back("size");
  // meta_name_.push_back("array_size");
  // meta_name_.push_back("num_fields");
  // meta_name_.push_back("ghosts_allocated");
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
  p | include_ghosts_;
}

//----------------------------------------------------------------------

void IoFieldData::meta_value
(int index,
 void ** buffer, std::string * name, int * type,
 int * mx, int * my, int * mz) throw()
{
}

//----------------------------------------------------------------------

void IoFieldData::field_array
(void ** buffer, std::string * name, int * type,
 int * pmx, int * pmy, int * pmz,
 int * pnx, int * pny, int * pnz) throw()
/// @param buffer [out] pointer to start of memory containing the requested data item
/// @param name   [out] identifier for the requested data item
/// @param type   [out] type and precision of the the requested data item
/// @param pmx,pmy,pmz [out] pointer to the array data item's allocated size
/// @param pnx,pny,pnz [out] pointer to the array data item's actual size
{
  FieldDescr * field_descr = cello::field_descr();
  
  if (buffer) (*buffer) = (void * ) 
		field_data_->values(field_descr,field_index_);
  if (name)   (*name) = 
		std::string("field_") +	field_descr->field_name(field_index_);
  int type_size = 0;
  if (type) {

    precision_type precision = field_descr->precision(field_index_);
    if (precision == precision_default) precision = default_precision;
    switch (precision) {
    case precision_single:
      (*type) = type_float;
      type_size = sizeof(float);
      break;
    case precision_double:
      (*type) = type_double;
      type_size = sizeof(double);
      break;
    case precision_quadruple:
      (*type) = type_quadruple;
      type_size = sizeof(long double);
      break;
    default:
      ERROR2 ("IoFieldData",
	      "Unsupported precision type %d for field %s",
	      precision, field_descr->field_name(field_index_).c_str());
    }
  }
  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;
  int cx,cy,cz;

  field_data_->dimensions(field_descr,field_index_,&mx,&my,&mz);
  field_data_->size(&nx,&ny,&nz);
  field_descr->ghost_depth(field_index_,&gx,&gy,&gz);
  field_descr->centering(field_index_,&cx,&cy,&cz);

  if (pmx) (*pmx) = mx;
  if (pmy) (*pmy) = my;
  if (pmz) (*pmz) = mz;

  if ((! include_ghosts_) && (field_data_->ghosts_allocated())) {

    // adjust buffer pointer to start of non-ghost values
    if (buffer) (*buffer) += type_size*(gx+mx*(gy+my*gz));
    if (pnx) (*pnx) = nx + cx;
    if (pny) (*pny) = ny + cy;
    if (pnz) (*pnz) = nz + cz;

  } else {

    if (pnx) (*pnx) = mx;
    if (pny) (*pny) = my;
    if (pnz) (*pnz) = mz;

  }
}

//----------------------------------------------------------------------

