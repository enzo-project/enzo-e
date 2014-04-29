// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed May 19 18:17:50 PDT 2010
/// @brief    Implementation of the FieldBlock class

#include "cello.hpp"
#include "field.hpp"

//----------------------------------------------------------------------

FieldBlock::FieldBlock 
(
 FieldDescr * field_descr,
 int nx, int ny, int nz 
 ) throw()
  : field_descr_(field_descr),
    array_(),
    offsets_(),
    ghosts_allocated_(true)
{
  if (nx != 0) {
    size_[0] = nx;
    size_[1] = ny;
    size_[2] = nz;
  } else {
    size_[0] = 0;
    size_[1] = 0;
    size_[2] = 0;
  }
}

//----------------------------------------------------------------------

FieldBlock::~FieldBlock() throw()
{  
  deallocate_array();
}

//----------------------------------------------------------------------

FieldBlock::FieldBlock ( const FieldBlock & field_block ) throw ()
{
  
  INCOMPLETE("FieldBlock::FieldBlock"); 
}

//----------------------------------------------------------------------

FieldBlock & FieldBlock::operator= ( const FieldBlock & field_block ) throw ()
{  
  INCOMPLETE("FieldBlock::operator=");
  return *this;
}

//----------------------------------------------------------------------

void FieldBlock::pup(PUP::er &p)
{
  TRACEPUP;

  bool up = p.isUnpacking();

  if (up) field_descr_ = new FieldDescr;
  p | *field_descr_;

  PUParray(p,size_,3);

  p | array_;
  p | offsets_;
  p | ghosts_allocated_;
}


//----------------------------------------------------------------------

void FieldBlock::size( int * nx, int * ny, int * nz ) const throw()
{
  if (nx) (*nx) = size_[0];
  if (ny) (*ny) = size_[1];
  if (nz) (*nz) = size_[2];
}

//----------------------------------------------------------------------

char * FieldBlock::field_values ( int id_field ) 
  throw (std::out_of_range)
{
  char * field = 0;
  if (unsigned(id_field) < offsets_.size()) {
    field = &array_[0] + offsets_[id_field];
  }
  ASSERT3 ("FieldBlock::field_values()",
	   "Trying to access invalid field id %d out of %d in FieldBlock %p",
	   id_field,offsets_.size(),this,
	   field != NULL);
  return field;
}

//----------------------------------------------------------------------

const char * FieldBlock::field_values ( int id_field ) const 
  throw (std::out_of_range)
{
  return (unsigned(id_field) < offsets_.size()) ? 
    &array_[0] + offsets_[id_field] : NULL;
}

//----------------------------------------------------------------------

const char * FieldBlock::field_unknowns ( int id_field ) const
  throw (std::out_of_range)
{
  return (const char *)
    ((FieldBlock *)this) -> field_unknowns(id_field);
}

//----------------------------------------------------------------------

char * FieldBlock::field_unknowns ( int id_field  )
  throw (std::out_of_range)
{
  char * field_unknowns = &array_[0] + offsets_[id_field];

  if (field_descr_ == NULL) return NULL;

  if ( ghosts_allocated() ) {

    int gx,gy,gz;
    field_descr_->ghosts(id_field,&gx,&gy,&gz);

    bool cx,cy,cz;
    field_descr_->centering(id_field,&cx,&cy,&cz);

    int nx,ny,nz;
    size(&nx,&ny,&nz);

    nx += 2*gx + (cx?0:1);
    ny += 2*gy + (cy?0:1);
    nz += 2*gz + (cz?0:1);

    precision_type precision = field_descr_->precision(id_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    field_unknowns += bytes_per_element * (gx + nx*(gy + ny*gz));
  } 

  return field_unknowns;
}

//----------------------------------------------------------------------

void FieldBlock::cell_width
(
 double xm, double xp, double * hx,
 double ym, double yp, double * hy,
 double zm, double zp, double * hz
 ) const throw ()
{
  if (hx) (*hx) = (xp-xm) / size_[0];
  if (hy) (*hy) = (yp-ym) / size_[1];
  if (hz) (*hz) = (zp-zm) / size_[2];
}

//----------------------------------------------------------------------

void FieldBlock::clear
(
 float              value,
 int                id_field_first,
 int                id_field_last) throw()
{
  if (field_descr_ == NULL) {
    WARNING ("FieldBlock::clear()",
	     "Trying to clear an unallocated FieldBlock");
    return;
  }

  if ( array_allocated() ) {
    if (id_field_first == -1) {
      id_field_first = 0;
      id_field_last  = field_descr_->field_count() - 1;
    } else if (id_field_last == -1) {
      id_field_last  = id_field_first;
    }
    for (int id_field = id_field_first;
	 id_field <= id_field_last;
	 id_field++) {
      int nx,ny,nz;
      field_size(id_field,&nx,&ny,&nz);
      precision_type precision = field_descr_->precision(id_field);
      char * array = &array_[0] + offsets_[id_field];
      switch (precision) {
      case precision_single:
	for (int i=0; i<nx*ny*nz; i++) {
	  ((float *) array)[i] = (float) value;
	}
	break;
      case precision_double:
	for (int i=0; i<nx*ny*nz; i++) {
	  ((double *) array)[i] = (double) value;
	}
	break;
      case precision_quadruple:
	for (int i=0; i<nx*ny*nz; i++) {
	  ((long double *) array)[i] = (long double) value;
	}
	break;
      default:
	char buffer[80];
	sprintf (buffer,"Clear called with unsupported precision %s" , 
		 cello::precision_name[precision]);
	ERROR("FieldBlock::clear", buffer);
      }

    }
  } else {
    ERROR("FieldBlock::clear",
	  "Called clear with unallocated arrays");
  }
}

//----------------------------------------------------------------------

void FieldBlock::allocate_array ( bool ghosts_allocated ) throw()
{

  // Error check size

  if (! (size_[0] > 0 &&
	 size_[1] > 0 &&
	 size_[2] > 0) ) {
    ERROR ("FieldBlock::allocate_array",
		   "Allocate called with zero field size");
  }

  // Warning check array already allocated
  if (array_allocated() ) {
    WARNING ("FieldBlock::allocate_array",
	     "Array already allocated: calling reallocate()");
    reallocate_array(ghosts_allocated);
    return;
  }

  ghosts_allocated_ = ghosts_allocated;

  int padding   = field_descr_->padding();
  int alignment = field_descr_->alignment();

  int array_size = 0;

  for (int id_field=0; id_field<field_descr_->field_count(); id_field++) {

    // Increment array_size, including padding and alignment adjustment

    int nx,ny,nz;       // not needed

    int size = field_size(id_field, &nx,&ny,&nz);

    array_size += adjust_padding_   (size,padding);
    array_size += adjust_alignment_ (size,alignment);

  }

  // Adjust for possible initial misalignment

  array_size += alignment - 1;

  // Allocate the array

  array_.resize(array_size);

  // Initialize field_begin

  int field_offset = align_padding_(alignment);

  offsets_.reserve(field_descr_->field_count());

  for (int id_field=0; id_field<field_descr_->field_count(); id_field++) {

    offsets_.push_back(field_offset);

    // Increment array_size, including padding and alignment adjustment

    int nx,ny,nz;       // not needed

    int size = field_size(id_field,&nx,&ny,&nz);

    field_offset += adjust_padding_  (size,padding);
    field_offset += adjust_alignment_(size,alignment);
  }

  // check if array_size is too big or too small

  if ( ! ( 0 <= (array_size - field_offset)
	   &&   (array_size - field_offset) < alignment)) {
    ERROR ("FieldBlock::allocate_array",
	   "Code error: array size was computed incorrectly");
  };

}

//----------------------------------------------------------------------

void FieldBlock::reallocate_array
(
 bool               ghosts_allocated
 ) throw()
{
  if (! array_allocated() ) {
    WARNING ("FieldBlock::reallocate_array",
	     "Array not allocated yet: calling allocate()");
    allocate_array(ghosts_allocated);
    return;
  }
  
  std::vector<int>  old_offsets;
  std::vector<char> old_array;

  old_array = array_;
  old_offsets = offsets_;

  array_.clear();
  offsets_.clear();

  ghosts_allocated_ = ghosts_allocated;

  allocate_array(ghosts_allocated_);

  restore_array_ (&old_array[0], old_offsets);
}

//----------------------------------------------------------------------

void FieldBlock::deallocate_array () throw()
{
  if ( array_allocated() ) {

    array_.clear();
    offsets_.clear();
  }
  // if (field_faces_ != 0) {
  //   delete field_faces_;
  //   field_faces_ = 0;
  // }

}

//----------------------------------------------------------------------

// void FieldBlock::allocate_ghosts(FieldDescr * field_descr) throw ()
// {
//   if (! ghosts_allocated() ) {

//   } else {
//     WARNING("FieldBlock::allocate_ghosts",
// 		    "Allocate called with ghosts already allocated");
//   }
// }

// //----------------------------------------------------------------------

// void FieldBlock::deallocate_ghosts(FieldDescr * field_descr) throw ()
// {
//   if ( ghosts_allocated() ) {

//     std::vector<int> old_offsets;
//     std::vector<char> old_array;

//     old_array = array_;

//     backup_array_ (field_descr,old_offsets);

//     ghosts_allocated_ = false;

//     allocate_array(field_descr);

//     restore_array_ (field_descr, &old_array[0], old_offsets);

//   } else {
//     WARNING("FieldBlock::deallocate_ghosts",
// 	    "Function called with ghosts not allocated");
//   }
// }

//======================================================================

int FieldBlock::adjust_padding_
(
 int size, 
 int padding) const throw ()
{
  return size + padding;
}

//----------------------------------------------------------------------

int FieldBlock::adjust_alignment_
(
 int size, 
 int alignment) const throw ()
{
  return (alignment - (size % alignment)) % alignment;
}

//----------------------------------------------------------------------

int FieldBlock::align_padding_ (int alignment) const throw()
{ 
  long unsigned start_long = reinterpret_cast<long unsigned>(&array_[0]);
  return ( alignment - (start_long % alignment) ) % alignment; 
}

//----------------------------------------------------------------------

int FieldBlock::field_size
(
 int                id_field,
 int              * nx,
 int              * ny,
 int              * nz
 ) const throw()
{

  if (field_descr_ == NULL) {
    WARNING("FieldBlock::field_size",
 	    "Calling field_size() on unallocated FieldBlock");
    if (nx) (*nx) = 0;
    if (ny) (*ny) = 0;
    if (nz) (*nz) = 0;
  }

  // Adjust memory usage due to ghosts if needed

  int  gx,gy,gz;
  if ( ghosts_allocated_ ) {
    field_descr_->ghosts(id_field,&gx,&gy,&gz);
  } else {
    gx = gy = gz = 0;
  }

  // Adjust memory usage due to field centering if needed

  bool cx,cy,cz;
  field_descr_->centering(id_field,&cx,&cy,&cz);

  // Compute array size

  if (nx) (*nx) = size_[0] + (1-cx) + 2*gx;
  if (ny) (*ny) = size_[1] + (1-cy) + 2*gy;
  if (nz) (*nz) = size_[2] + (1-cz) + 2*gz;

  // Return array size in bytes

  precision_type precision = field_descr_->precision(id_field);
  int bytes_per_element = cello::sizeof_precision (precision);
  if (nz) {
    return (*nx) * (*ny) * (*nz) * bytes_per_element;
  } else if (ny) {
    return (*nx) * (*ny) * bytes_per_element;
  } else {
    return (*nx) * bytes_per_element;
  }
}

//----------------------------------------------------------------------

void FieldBlock::mul (int id_1, int id_2)
{
  if (field_descr_ == NULL) {
    WARNING("FieldBlock::mul()",
 	    "Calling operation on unallocated FieldBlock");
    return;
  }

  precision_type p1 = field_descr_->precision(id_1);
  precision_type p2 = field_descr_->precision(id_2);
  ASSERT ("FieldBlock::mul",
	  "FieldBlock precisions must be constant",
	  p1 == p2);

  char * field_1 = field_unknowns(id_1);
  char * field_2 = field_unknowns(id_2);
  
  switch (p1) {
  case precision_single:
    mul_((float *) field_1, (float *) field_2);
    break;
  case precision_double:
    mul_((double *) field_1,(double *) field_2);
    break;
  default:
    break;
  }

}

//----------------------------------------------------------------------

template<class T>
void FieldBlock::mul_ (T * field_1, T * field_2)
{
  int nx,ny,nz;
  size (&nx,&ny,&nz);
  int gx,gy,gz;
  field_descr_->ghosts (0,&gx,&gy,&gz);
  int ndx = (nx>1) ? nx+2*gx : 1;
  int ndy = (ny>1) ? ny+2*gy : 1;

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + ndx*(iy + ndy*iz);
	field_1[i] *= field_2[i];
      }
    }
  }
}

//----------------------------------------------------------------------

void FieldBlock::div (int id_1, int id_2)
{
  precision_type p1 = field_descr_->precision(id_1);
  precision_type p2 = field_descr_->precision(id_2);
  ASSERT ("FieldBlock::div",
	  "FieldBlock precisions must be constant",
	  p1 == p2);

  char * field_1 = field_unknowns(id_1);
  char * field_2 = field_unknowns(id_2);
  
  switch (p1) {
  case precision_single:
    div_((float *) field_1, (float *) field_2);
    break;
  case precision_double:
    div_((double *) field_1,(double *) field_2);
    break;
  default:
    break;
  }

}


//----------------------------------------------------------------------

void FieldBlock::scale (int id, double value)
{
  precision_type p = field_descr_->precision(id);

  char * field = field_unknowns(id);
  
  switch (p) {
  case precision_single:
    scale_((float *) field, value);
    break;
  case precision_double:
    scale_((double *) field,value);
    break;
  default:
    break;
  }

}


//----------------------------------------------------------------------

template<class T>
void FieldBlock::div_ (T * field_1, T * field_2)
{
  int nx,ny,nz;
  size (&nx,&ny,&nz);
  int gx,gy,gz;
  field_descr_->ghosts (0,&gx,&gy,&gz);
  int ndx = (nx>1) ? nx+2*gx : 1;
  int ndy = (ny>1) ? ny+2*gy : 1;

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + ndx*(iy + ndy*iz);
	field_1[i] /= field_2[i];
      }
    }
  }
}

//----------------------------------------------------------------------

template<class T>
void FieldBlock::scale_ (T * field, double value)
{
  int nx,ny,nz;
  size (&nx,&ny,&nz);
  int gx,gy,gz;
  field_descr_->ghosts (0,&gx,&gy,&gz);
  int ndx = (nx>1) ? nx+2*gx : 1;
  int ndy = (ny>1) ? ny+2*gy : 1;

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + ndx*(iy + ndy*iz);
	field[i] *= value;
      }
    }
  }
}

//----------------------------------------------------------------------

void FieldBlock::print
(const char * message,
 // double lower[3],
 // double upper[3],
 bool use_file) const throw()
{

#ifndef CELLO_DEBUG
  return;
#else

  int ip=0;

   ip=CkMyPe();

   char filename [80];
   sprintf (filename,"%s-%d.debug",message,ip);
   printf ("DEBUG message = %s\n",message);
   printf ("DEBUG filename = %s\n",filename);

   FILE * fp = fopen (filename,"w");

   ASSERT("FieldBlock::print",
	  "FieldBlocks not allocated",
	  array_allocated());

   int field_count = field_descr_->field_count();
   for (int index_field=0; index_field<field_count; index_field++) {

     // WARNING: not copying string works on some compilers but not others
     const char * field_name = strdup(field_descr_->field_name(index_field).c_str());

     int nxd,nyd,nzd;
     field_size(index_field,&nxd,&nyd,&nzd);
     int gx,gy,gz;
     field_descr_->ghosts(index_field,&gx,&gy,&gz);

     int ixm,iym,izm;
     int ixp,iyp,izp;

     // Exclude ghost zones

     // ixm = gx;
     // iym = gy;
     // izm = gz;

     // ixp = nxd - gx;
     // iyp = nyd - gy;
     // izp = nzd - gz;

     // Include ghost zones

     ixm = 0;
     iym = 0;
     izm = 0;

     ixp = nxd;
     iyp = nyd;
     izp = nzd;

     int nx,ny,nz;

     nx = (ixp-ixm);
     ny = (iyp-iym);
     nz = (izp-izm);

     //     double hx,hy,hz;

     // hx = (upper[0]-lower[0])/(nxd-2*gx);
     // hy = (upper[1]-lower[1])/(nyd-2*gy);
     // hz = (upper[2]-lower[2])/(nzd-2*gz);

     const char * array_offset = &array_[0]+offsets_[index_field];
     switch (field_descr->precision(index_field)) {
     case precision_single:
       print_((const float * ) array_offset,
	      field_name, message, // lower,
	      fp,
	      ixm,iym,izm,
	      ixp,iyp,izp,
	      nx, ny, nz,
	      gx, gy ,gz,
	      //	      hx, hy ,hz,
	      nxd,nyd);
       break;
     case precision_double:
       print_((const double * ) array_offset, 
	      field_name, message, // lower,
	      fp,
	      ixm,iym,izm,
	      ixp,iyp,izp,
	      nx, ny, nz,
	      gx, gy ,gz,
	      //	      hx, hy ,hz,
	      nxd,nyd);
       break;
     case precision_quadruple:
       print_((const long double * ) array_offset, 
	      field_name, message, // lower,
	      fp,
	      ixm,iym,izm,
	      ixp,iyp,izp,
	      nx, ny, nz,
	      gx, gy ,gz,
	      //	      hx, hy ,hz,
	      nxd,nyd);
       break;
     default:
       ERROR("FieldBlock::print", "Unsupported precision");
     }

     free ((void *)field_name);
   }

#endif /* ifndef CELLO_DEBUG */
}

template <class T>
void FieldBlock::print_
(const T * field,
 const char * field_name,
 const char * message,
 // double lower [3],
 FILE * fp,
 int ixm,int iym,int izm,
 int ixp,int iyp,int izp,
 int nx, int ny, int nz,
 int gx, int gy ,int gz,
 // double hx, double hy ,double hz,
 int nxd,int nyd) const
{

  T min = std::numeric_limits<T>::max();
  T max = std::numeric_limits<T>::min();
  double sum = 0.0;
  for (int iz=izm; iz<izp; iz++) {
    for (int iy=iym; iy<iyp; iy++) {
      for (int ix=ixm; ix<ixp; ix++) {
	int i = ix + nxd*(iy + nyd*iz);
	min = MIN(min,field[i]);
	max = MAX(max,field[i]);
	sum += field[i];
#ifdef CELLO_DEBUG_VERBOSE
	// double x = hx*(ix-gx) + lower[axis_x];
	// double y = hy*(iy-gy) + lower[axis_y];
	// double z = hz*(iz-gz) + lower[axis_z];
	if (isnan(field[i])) {
	  fprintf(fp,"DEBUG: %s %s  %2d %2d %2d NAN\n",
		  message,field_name,ix,iy,iz);
	} else {
	  fprintf(fp,"DEBUG: %s %s  %2d %2d %2d %f\n",
		  message,field_name,ix,iy,iz,field[i]);
	}
#endif
      }
    }
  }
  double avg = sum / (nx*ny*nz);
  fprintf
    (fp,"%s [%s] %18.14g %18.14g %18.14g\n",
     message ? message : "", field_name, min,avg,max);

}

//----------------------------------------------------------------------

void FieldBlock::restore_array_ 
( const char * array_from,
  std::vector<int> & offsets_from) throw (std::out_of_range)
{

  // copy values
  for (int id_field=0; 
       id_field < field_descr_->field_count();
       id_field++) {

    // get "to" field size

    int nx2,ny2,nz2;
    field_size(id_field, &nx2,&ny2,&nz2);

    // get "from" field size

    ghosts_allocated_ = ! ghosts_allocated_;

    int nx1,ny1,nz1;
    field_size(id_field, &nx1,&ny1,&nz1);

    ghosts_allocated_ = ! ghosts_allocated_;

    // determine offsets to unknowns if ghosts allocated

    int offset1 = (nx1-nx2)/2 + nx1* ( (ny1-ny2)/2 + ny1 * (nz1-nz2)/2 );
    offset1 = MAX (offset1, 0);

    int offset2 = (nx2-nx1)/2 + nx2* ( (ny2-ny1)/2 + ny2 * (nz2-nz1)/2 );
    offset2 = MAX (offset2, 0);

    // determine unknowns size

    int nx = MIN(nx1,nx2);
    int ny = MIN(ny1,ny2);
    int nz = MIN(nz1,nz2);

    // adjust for precision

    precision_type precision = field_descr_->precision(id_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    offset1 *= bytes_per_element;
    offset2 *= bytes_per_element;

    // determine array start

    const char * array1 = array_from + offsets_from.at(id_field) + offset1;
    char       * array2 = &array_[0] + offsets_.at(id_field)     + offset2;

    // copy values (use memcopy?)

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  for (int ip=0; ip<bytes_per_element; ip++) {
	    int i1 = ip + bytes_per_element*(ix + nx1*(iy + ny1*iz));
	    int i2 = ip + bytes_per_element*(ix + nx2*(iy + ny2*iz));
	    array2[i2] = array1[i1]; // XXX array1 = garbage
	  }
	}
      }
    }
  }

  offsets_from.clear();
}
