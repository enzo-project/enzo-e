// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed May 19 18:17:50 PDT 2010
/// @brief    Implementation of the FieldData class

#include "cello.hpp"
#include "field.hpp"

#include "charm_simulation.hpp"

//----------------------------------------------------------------------

FieldData::FieldData 
(
 FieldDescr * field_descr,
 int nx, int ny, int nz 
 ) throw()
  : field_descr_(field_descr),
    array_permanent_(),
    array_temporary_(),
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

FieldData::~FieldData() throw()
{  
  deallocate_permanent();
}

//----------------------------------------------------------------------

FieldData::FieldData ( const FieldData & field_data ) throw ()
{
  
  INCOMPLETE("FieldData::FieldData"); 
}

//----------------------------------------------------------------------

FieldData & FieldData::operator= ( const FieldData & field_data ) throw ()
{  
  INCOMPLETE("FieldData::operator=");
  return *this;
}

//----------------------------------------------------------------------

void FieldData::pup(PUP::er &p)
{
  TRACEPUP;

  bool up = p.isUnpacking();

  if (up) field_descr_ = new FieldDescr;
  p | *field_descr_;

  PUParray(p,size_,3);

  p | array_permanent_;
  //  p | array_temporary_;
  static bool warn = false;
  if (! warn) {
    WARNING("FieldData::pup()",
	    "Skipping array_temporary_");
    warn = true;
  }
  p | offsets_;
  p | ghosts_allocated_;
}


//----------------------------------------------------------------------

void FieldData::dimensions( int id_field, int * mx, int * my, int * mz ) const throw()
{
  int nx,ny,nz;
  int gx,gy,gz;
  int cx,cy,cz;
  size      (&nx,&ny,&nz);
  ghosts    (id_field,&gx,&gy,&gz);
  centering (id_field,&cx,&cy,&cz);

  if (mx) (*mx) = (nx > 1) ? (nx + 2*gx + cx) : 1;
  if (my) (*my) = (ny > 1) ? (ny + 2*gy + cy) : 1;
  if (mz) (*mz) = (nz > 1) ? (nz + 2*gz + cz) : 1;
}

//----------------------------------------------------------------------

void FieldData::size( int * nx, int * ny, int * nz ) const throw()
{
  if (nx) (*nx) = size_[0];
  if (ny) (*ny) = size_[1];
  if (nz) (*nz) = size_[2];
}

//----------------------------------------------------------------------

const char * FieldData::values ( int id_field ) const 
  throw (std::out_of_range)
{
  return (const char *)
    ((FieldData *)this) -> values(id_field);
}

//----------------------------------------------------------------------

char * FieldData::values ( int id_field ) 
  throw (std::out_of_range)
{
  char * values = 0;
  if (is_permanent(id_field)) {

    // permanent field

    if (0 <= id_field && id_field < field_count()) {
      values = &array_permanent_[0] + offsets_[id_field];
    }
  } else {

    // temporary field

    int id_temporary = id_field - num_permanent();

    if (0 <= id_temporary && id_temporary < array_temporary_.size()) {
      values = array_temporary_[id_temporary];
    }
  }
  return values;
}

//----------------------------------------------------------------------

const char * FieldData::unknowns ( int id_field ) const
  throw (std::out_of_range)
{
  return (const char *)
    ((FieldData *)this) -> unknowns(id_field);
}

//----------------------------------------------------------------------

char * FieldData::unknowns ( int id_field  )
  throw (std::out_of_range)
{
  char * unknowns = values(id_field);

  if ( ghosts_allocated() && unknowns ) {

    int gx,gy,gz;
    int mx,my;

    ghosts    (id_field,&gx,&gy,&gz);
    dimensions(id_field,&mx,&my);

    precision_type precision = this->precision(id_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    unknowns += bytes_per_element * (gx + mx*(gy + my*gz));
  } 
  return unknowns;
}

//----------------------------------------------------------------------

void FieldData::cell_width
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

void FieldData::clear
(
 float              value,
 int                id_field_first,
 int                id_field_last) throw()
{
  if (field_descr_ == NULL) {
    WARNING ("FieldData::clear()",
	     "Trying to clear an unallocated FieldData");
    return;
  }

  if ( permanent_allocated() ) {
    if (id_field_first == -1) {
      id_field_first = 0;
      id_field_last  = field_count() - 1;
    } else if (id_field_last == -1) {
      id_field_last  = id_field_first;
    }
    for (int id_field = id_field_first;
	 id_field <= id_field_last;
	 id_field++) {
      int nx,ny,nz;
      field_size(id_field,&nx,&ny,&nz);
      precision_type precision = this->precision(id_field);
      char * array = &array_permanent_[0] + offsets_[id_field];
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
	ERROR("FieldData::clear", buffer);
      }

    }
  } else {
    ERROR("FieldData::clear",
	  "Called clear with unallocated arrays");
  }
}

//----------------------------------------------------------------------

void FieldData::allocate_permanent ( bool ghosts_allocated ) throw()
{

  // Error check size

  if (! (size_[0] > 0 &&
	 size_[1] > 0 &&
	 size_[2] > 0) ) {
    ERROR ("FieldData::allocate_permanent",
		   "Allocate called with zero field size");
  }

  // Warning check array already allocated
  if (permanent_allocated() ) {
    WARNING ("FieldData::allocate_permanent",
	     "Array already allocated: calling reallocate()");
    reallocate_permanent(ghosts_allocated);
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

  array_permanent_.resize(array_size);

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
    ERROR ("FieldData::allocate_permanent",
	   "Code error: array size was computed incorrectly");
  };

}

//----------------------------------------------------------------------

void FieldData::allocate_temporary (int id_field) throw (std::out_of_range)

{
  int index_field = id_field - num_permanent();

  if (! (index_field < array_temporary_.size())) {
    array_temporary_.resize(index_field+1, 0);
  }
    
  if (array_temporary_[index_field] == 0) {
    int mx,my,mz;
    dimensions(id_field,&mx,&my,&mz);
    int m = mx*my*mz;
    precision_type precision = this->precision(id_field);
    if (precision == precision_single)    
      array_temporary_[index_field] = (char*) new float [m];
    if (precision == precision_double)    
      array_temporary_[index_field] = (char*) new double [m];
    if (precision == precision_quadruple) 
      array_temporary_[index_field] = (char*) new long double [m];
  } else {
    WARNING("FieldData::allocate_temporary",
	    "Calling allocate_temporary() on already-allocated Field");
  }
}

//----------------------------------------------------------------------

void FieldData::deallocate_temporary (int id_field) throw(std::out_of_range)
{
  int index_field = id_field - num_permanent();

  if (! (index_field < array_temporary_.size())) {
    array_temporary_.resize(index_field+1, 0);
  }
  if (array_temporary_[index_field] != 0) {
    precision_type precision = this->precision(id_field);
    if (precision == precision_single)    
      delete [] (float *)       array_temporary_[index_field];
    if (precision == precision_double)    
      delete [] (double *)      array_temporary_[index_field];
    if (precision == precision_quadruple) 
      delete [] (long double *) array_temporary_[index_field];
  }
  array_temporary_[index_field] = 0;

}
//----------------------------------------------------------------------

void FieldData::reallocate_permanent
(
 bool               ghosts_allocated
 ) throw()
{
  if (! permanent_allocated() ) {
    WARNING ("FieldData::reallocate_permanent",
	     "Array not allocated yet: calling allocate()");
    allocate_permanent(ghosts_allocated);
    return;
  }
  
  std::vector<int>  old_offsets;
  std::vector<char> old_array;

  old_array = array_permanent_;
  old_offsets = offsets_;

  array_permanent_.clear();
  offsets_.clear();

  ghosts_allocated_ = ghosts_allocated;

  allocate_permanent(ghosts_allocated_);

  restore_permanent_ (&old_array[0], old_offsets);
}

//----------------------------------------------------------------------

void FieldData::deallocate_permanent () throw()
{
  if ( permanent_allocated() ) {

    array_permanent_.clear();
    offsets_.clear();
  }
  // if (field_faces_ != 0) {
  //   delete field_faces_;
  //   field_faces_ = 0;
  // }

}

//----------------------------------------------------------------------

// void FieldData::allocate_ghosts(FieldDescr * field_descr) throw ()
// {
//   if (! ghosts_allocated() ) {

//   } else {
//     WARNING("FieldData::allocate_ghosts",
// 		    "Allocate called with ghosts already allocated");
//   }
// }

// //----------------------------------------------------------------------

// void FieldData::deallocate_ghosts(FieldDescr * field_descr) throw ()
// {
//   if ( ghosts_allocated() ) {

//     std::vector<int> old_offsets;
//     std::vector<char> old_array;

//     old_array = array_permanent_;

//     backup_array_ (field_descr,old_offsets);

//     ghosts_allocated_ = false;

//     allocate_array(field_descr);

//     restore_array_ (field_descr, &old_array[0], old_offsets);

//   } else {
//     WARNING("FieldData::deallocate_ghosts",
// 	    "Function called with ghosts not allocated");
//   }
// }

//======================================================================

int FieldData::adjust_padding_
(
 int size, 
 int padding) const throw ()
{
  return size + padding;
}

//----------------------------------------------------------------------

int FieldData::adjust_alignment_
(
 int size, 
 int alignment) const throw ()
{
  return (alignment - (size % alignment)) % alignment;
}

//----------------------------------------------------------------------

int FieldData::align_padding_ (int alignment) const throw()
{ 
  long unsigned start_long = reinterpret_cast<long unsigned>(&array_permanent_[0]);
  return ( alignment - (start_long % alignment) ) % alignment; 
}

//----------------------------------------------------------------------

int FieldData::field_size
(
 int                id_field,
 int              * nx,
 int              * ny,
 int              * nz
 ) const throw()
{

  if (field_descr_ == NULL) {
    WARNING("FieldData::field_size",
 	    "Calling field_size() on unallocated FieldData");
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

  int cx,cy,cz;
  field_descr_->centering(id_field,&cx,&cy,&cz);

  // Compute array size

  if (nx) (*nx) = size_[0] + 2*gx + cx;
  if (ny) (*ny) = size_[1] + 2*gy + cy;
  if (nz) (*nz) = size_[2] + 2*gz + cz;

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

void FieldData::print
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

   ASSERT("FieldData::print",
	  "FieldData not allocated",
	  permanent_allocated());

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

     const char * array_offset = &array_permanent_[0]+offsets_[index_field];
     switch (field_descr_->precision(index_field)) {
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
       ERROR("FieldData::print", "Unsupported precision");
     }

     free ((void *)field_name);
   }

#endif /* ifndef CELLO_DEBUG */
}

template <class T>
void FieldData::print_
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

void FieldData::restore_permanent_ 
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

    const char * array1 = offset1 + offsets_from.at(id_field) + array_from;
    char       * array2 = offset2 + offsets_.at    (id_field) + &array_permanent_[0];

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
