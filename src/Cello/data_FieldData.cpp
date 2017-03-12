// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed May 19 18:17:50 PDT 2010
/// @brief    Implementation of the FieldData class

#include "cello.hpp"
#include "data.hpp"

#include "charm_simulation.hpp"

//----------------------------------------------------------------------

FieldData::FieldData 
(
 const FieldDescr * field_descr,
 int nx, int ny, int nz 
 ) throw()
  : array_permanent_(),
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

void FieldData::pup(PUP::er &p)
{
  TRACEPUP;

  PUParray(p,size_,3);

  p | array_permanent_;
  //  p | array_temporary_;
  static bool warn[CONFIG_NODE_SIZE] = {false};
  const int in = cello::index_static();
  if (! warn[in]) {
    WARNING("FieldData::pup()",
  	    "Skipping array_temporary_");
    warn[in] = true;
  }
  p | offsets_;
  p | ghosts_allocated_;
}


//----------------------------------------------------------------------

void FieldData::dimensions(const FieldDescr * field_descr,
			   int id_field, int * mx, int * my, int * mz ) const throw()
{
  int nx,ny,nz;
  int gx,gy,gz;
  int cx,cy,cz;
  size      (&nx,&ny,&nz);
  field_descr->ghost_depth (id_field,&gx,&gy,&gz);
  field_descr->centering (id_field,&cx,&cy,&cz);

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

const char * FieldData::values ( const FieldDescr * field_descr,
				 int id_field ) const 
  throw ()
{
  return (const char *)
    ((FieldData *)this) -> values(field_descr,id_field);
}

//----------------------------------------------------------------------

char * FieldData::values (const FieldDescr * field_descr,
			  int id_field ) 
  throw ()
{
  if (id_field == -1) return NULL;

  char * values = 0;

  if (field_descr->is_permanent(id_field)) {

    // permanent field
    const int num_fields = field_descr->field_count();
    if (0 <= id_field && id_field < num_fields) {
      values = &array_permanent_[0] + offsets_[id_field];
    }

  } else {

    // temporary field

    int id_temporary = id_field - field_descr->num_permanent();

    if (0 <= id_temporary && id_temporary < int(array_temporary_.size())) {
      values = array_temporary_[id_temporary];
    }
  }
  return values;
}

//----------------------------------------------------------------------

const char * FieldData::unknowns ( const FieldDescr * field_descr, 
				   int id_field ) const
  throw ()
{
  return (const char *)
    ((FieldData *)this) -> unknowns(field_descr,id_field);
}

//----------------------------------------------------------------------

char * FieldData::unknowns (const FieldDescr * field_descr,
			    int id_field  )
  throw ()
{
  char * unknowns = values(field_descr,id_field);

  if ( ghosts_allocated() && unknowns ) {

    int gx,gy,gz;
    int mx,my;

    field_descr->ghost_depth    (id_field,&gx,&gy,&gz);
    dimensions(field_descr,id_field,&mx,&my);

    precision_type precision = field_descr->precision(id_field);
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
 const FieldDescr * field_descr,
 float              value,
 int                id_field_first,
 int                id_field_last) throw()
{
  if ( permanent_allocated() ) {
    if (id_field_first == -1) {
      id_field_first = 0;
      id_field_last  = field_descr->field_count() - 1;
    } else if (id_field_last == -1) {
      id_field_last  = id_field_first;
    }
    for (int id_field = id_field_first;
	 id_field <= id_field_last;
	 id_field++) {
      int nx,ny,nz;
      field_size(field_descr,id_field,&nx,&ny,&nz);
      precision_type precision = field_descr->precision(id_field);
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

void FieldData::allocate_permanent 
(const FieldDescr * field_descr,
 bool ghosts_allocated ) throw()
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
    reallocate_permanent(field_descr,ghosts_allocated);
    return;
  }

  ghosts_allocated_ = ghosts_allocated;

  int padding   = field_descr->padding();
  int alignment = field_descr->alignment();

  int array_size = 0;

  for (int id_field=0; id_field<field_descr->field_count(); id_field++) {

    // Increment array_size, including padding and alignment adjustment

    int nx,ny,nz;       // not needed

    int size = field_size(field_descr,id_field, &nx,&ny,&nz);

    array_size += adjust_padding_   (size,padding);
    array_size += adjust_alignment_ (size,alignment);

  }

  // Adjust for possible initial misalignment

  array_size += alignment - 1;

  // Allocate the array

  array_permanent_.resize(array_size);

  // Initialize field_begin

  int field_offset = align_padding_(alignment);

  offsets_.reserve(field_descr->field_count());

  for (int id_field=0; id_field<field_descr->field_count(); id_field++) {

    offsets_.push_back(field_offset);

    // Increment array_size, including padding and alignment adjustment

    int nx,ny,nz;       // not needed

    int size = field_size(field_descr,id_field,&nx,&ny,&nz);

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

void FieldData::allocate_temporary (const FieldDescr * field_descr,
				    int id_field) throw ()

{
  int index_field = id_field - field_descr->num_permanent();

  if (! (index_field < int(array_temporary_.size()))) {
    array_temporary_.resize(index_field+1, 0);
  }
    
  if (array_temporary_[index_field] == 0) {
    int mx,my,mz;
    dimensions(field_descr,id_field,&mx,&my,&mz);
    int m = mx*my*mz;
    precision_type precision = field_descr->precision(id_field);
    if (precision == precision_single) {
      array_temporary_[index_field] = (char*) new float [m];
    } else if (precision == precision_double) {
      array_temporary_[index_field] = (char*) new double [m];
    } else if (precision == precision_quadruple) {
      array_temporary_[index_field] = (char*) new long double [m];
    } else {
      WARNING("FieldData::allocate_temporary",
	      "Calling allocate_temporary() on already-allocated Field");
    }
  }
}

//----------------------------------------------------------------------

void FieldData::deallocate_temporary (const FieldDescr * field_descr,
				      int id_field) throw()
{
  int index_field = id_field - field_descr->num_permanent();

  if (! (index_field < int(array_temporary_.size()))) {
    array_temporary_.resize(index_field+1, 0);
  }
  if (array_temporary_[index_field] != 0) {
    precision_type precision = field_descr->precision(id_field);
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
 const FieldDescr * field_descr,
 bool               ghosts_allocated
 ) throw()
{
  if (! permanent_allocated() ) {
    WARNING ("FieldData::reallocate_permanent",
	     "Array not allocated yet: calling allocate()");
    allocate_permanent(field_descr,ghosts_allocated);
    return;
  }
  
  std::vector<int>  old_offsets;
  std::vector<char> old_array;

  old_array = array_permanent_;
  old_offsets = offsets_;

  array_permanent_.clear();
  offsets_.clear();

  ghosts_allocated_ = ghosts_allocated;

  allocate_permanent(field_descr,ghosts_allocated_);

  restore_permanent_ (field_descr,&old_array[0], old_offsets);
}

//----------------------------------------------------------------------

void FieldData::deallocate_permanent () throw()
{
  if ( permanent_allocated() ) {

    array_permanent_.clear();
    offsets_.clear();
  }
}

//----------------------------------------------------------------------

int FieldData::field_size
(
 const FieldDescr * field_descr,
 int                id_field,
 int              * nx,
 int              * ny,
 int              * nz
 ) const throw()
{
  // Adjust memory usage due to ghosts if needed

  int  gx,gy,gz;
  if ( ghosts_allocated_ ) {
    field_descr->ghost_depth(id_field,&gx,&gy,&gz);
  } else {
    gx = gy = gz = 0;
  }

  // Adjust memory usage due to field centering if needed

  int cx,cy,cz;
  field_descr->centering(id_field,&cx,&cy,&cz);

  // Compute array size

  if (nx) (*nx) = size_[0] + 2*gx + cx;
  if (ny) (*ny) = size_[1] + 2*gy + cy;
  if (nz) (*nz) = size_[2] + 2*gz + cz;

  // Return array size in bytes

  precision_type precision = field_descr->precision(id_field);
  int bytes_per_element = cello::sizeof_precision (precision);

  int bytes_total = bytes_per_element;

  if (nx) bytes_total *= (*nx);
  if (ny) bytes_total *= (*ny);
  if (nz) bytes_total *= (*nz);

  return bytes_total;
}

//----------------------------------------------------------------------

void FieldData::print
(
 const FieldDescr * field_descr,
 const char * message,
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

   int field_count = field_descr->field_count();
   for (int index_field=0; index_field<field_count; index_field++) {

     // WARNING: not copying string works on some compilers but not others
     const char * field_name = strdup(field_descr->field_name(index_field).c_str());

     int nxd,nyd,nzd;
     field_size(field_descr,index_field,&nxd,&nyd,&nzd);
     int gx,gy,gz;
     field_descr->ghost_depth(index_field,&gx,&gy,&gz);

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
       ERROR("FieldData::print", "Unsupported precision");
     }

     free ((void *)field_name);
   }

#endif /* ifndef CELLO_DEBUG */
}

//----------------------------------------------------------------------

// void FieldData::copy
// (const FieldDescr * field_descr,
//  int id, int is,
//  bool ghosts = true ) throw()
// {
// }

// //----------------------------------------------------------------------

// void FieldData::axpy
// (const FieldDescr * field_descr, int iz, double a, int ix, int iy,
//  bool ghosts = true ) throw()
// {
// }

//----------------------------------------------------------------------

// template<class T>
//   void FieldData::axpy_ (T* Z, double a, const T* X, const T* Y) const throw()
// {
//   for (int iz=0; iz<mz_; iz++) 
//     for (int iy=0; iy<my_; iy++) 
//       for (int ix=0; ix<mx_; ix++) {
// 	int i = ix + mx_*(iy + my_*iz);
// 	Z[i] = a * X[i] + Y[i];
//       }
// }

//----------------------------------------------------------------------

void FieldData::axpy
(const FieldDescr * field_descr,
 int iz, long double a, int ix, int iy, bool ghosts) throw()
{
  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;

  size                        (&nx, &ny, &nz);
  dimensions  (field_descr,ix, &mx, &my, &mz);
  field_descr->ghost_depth(ix, &gx, &gy, &gz);

  void * x = values(field_descr,ix);
  void * y = values(field_descr,iy);
  void * z = values(field_descr,iz);

  switch (field_descr->precision(ix)) {
  case precision_single:
    axpy_ ((float*)z,a,(float*)x,(float*)y,ghosts,
	  mx,my,mz,nx,ny,nz,gx,gy,gz);
    break;
  case precision_double:
    axpy_ ((double*)z,a,(double*)x,(double*)y,ghosts,
	   mx,my,mz,nx,ny,nz,gx,gy,gz);
    break;
  case precision_quadruple:
    axpy_ ((long double*)z,a,(long double*)x,(long double*)y,ghosts,
	   mx,my,mz,nx,ny,nz,gx,gy,gz);
    break;
  default:
    ERROR2("FieldData::axpy()",
	   "Unknown precision %d for field id %d",
	   field_descr->precision(ix),ix);
    break;
  }

}

template<class T>
void FieldData::axpy_
(T * Z, long double a, const T * X, const T * Y, bool ghosts,
 int mx, int my, int mz,
 int nx, int ny, int nz,
 int gx, int gy, int gz) const throw()
{
  if (ghosts) {
    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
	for (int ix=0; ix<mx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  Z[i] = a * X[i] + Y[i];
	}
      }
    }
  } else {
    int i0 = gx + mx*(gy + my*gz);
    for (int iz=0; iz<nx; iz++) {
      for (int iy=0; iy<ny; iy++) { 
	for (int ix=0; ix<mx; ix++) {
	  int i = i0 + ix + mx*(iy + my*iz);
	  Z[i] = a * X[i] + Y[i];
	}
      }
    }
  }
}

//----------------------------------------------------------------------

double FieldData::dot (const FieldDescr * field_descr, int ix, int iy) throw()
{
  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;

  size                        (&nx, &ny, &nz);
  dimensions  (field_descr,ix, &mx, &my, &mz);
  field_descr->ghost_depth(ix, &gx, &gy, &gz);

  void * x = values(field_descr,ix);
  void * y = values(field_descr,iy);

  switch (field_descr->precision(ix)) {
  case precision_single:
    return dot_ ((float*)x,(float*)y,mx,my,mz,nx,ny,nz,gx,gy,gz);
    break;
  case precision_double:
    return dot_ ((double*)x,(double*)y,mx,my,mz,nx,ny,nz,gx,gy,gz);
    break;
  case precision_quadruple:
    return dot_ ((long double*)x,(long double*)y,mx,my,mz,nx,ny,nz,gx,gy,gz);
    break;
  default:
    ERROR2("FieldData::dot()",
	   "Unknown precision %d for field id %d",
	   field_descr->precision(ix),ix);
    return 0.0;
    break;
  }
}

template<class T>
long double FieldData::dot_(const T* X, const T* Y, int mx, int my, int mz, int nx, int ny, int nz, int gx, int gy, int gz) const throw()
{
  const int i0 = gx + mx*(gy + my*gz);
  long double value = 0.0;
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = i0 + (ix + mx*(iy + my*iz));
	value += X[i]*Y[i];
      }
    }
  }
  return value;
}

//----------------------------------------------------------------------

void FieldData::scale
(const FieldDescr * field_descr,
 int iy, long double a, int ix, bool ghosts) throw()
{
  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;

  size                        (&nx, &ny, &nz);
  dimensions  (field_descr,ix, &mx, &my, &mz);
  field_descr->ghost_depth(ix, &gx, &gy, &gz);

  void * x = values(field_descr,ix);
  void * y = values(field_descr,iy);

  switch (field_descr->precision(ix)) {
  case precision_single:
    scale_ ((float*)y,a,(float*)x,ghosts,
	  mx,my,mz,nx,ny,nz,gx,gy,gz);
    break;
  case precision_double:
    scale_ ((double*)y,a,(double*)x,ghosts,
	   mx,my,mz,nx,ny,nz,gx,gy,gz);
    break;
  case precision_quadruple:
    scale_ ((long double*)y,a,(long double*)x,ghosts,
	   mx,my,mz,nx,ny,nz,gx,gy,gz);
    break;
  default:
    ERROR2("FieldData::scale()",
	   "Unknown precision %d for field id %d",
	   field_descr->precision(ix),ix);
    break;
  }

}

template<class T>
void FieldData::scale_
(T * Y, long double a, T * X, bool ghosts,
 int mx, int my, int mz,
 int nx, int ny, int nz,
 int gx, int gy, int gz) const throw()
{
  if (ghosts) {
    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
	for (int ix=0; ix<mx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  Y[i] = a*X[i];
	}
      }
    }
  } else {
    int i0 = gx + mx*(gy + my*gz);
    for (int iz=0; iz<nx; iz++) {
      for (int iy=0; iy<ny; iy++) { 
	for (int ix=0; ix<mx; ix++) {
	  int i = i0 + ix + mx*(iy + my*iz);
	  Y[i] = a*X[i];
	}
      }
    }
  }
}

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
  T max = - std::numeric_limits<T>::max();
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
(
 const FieldDescr * field_descr,
 const char * array_from,
  std::vector<int> & offsets_from) throw ()
{

  // copy values
  for (int id_field=0; 
       id_field < field_descr->field_count();
       id_field++) {

    // get "to" field size

    int nx2,ny2,nz2;
    field_size(field_descr,id_field, &nx2,&ny2,&nz2);

    // get "from" field size

    ghosts_allocated_ = ! ghosts_allocated_;

    int nx1,ny1,nz1;
    field_size(field_descr,id_field, &nx1,&ny1,&nz1);

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

    precision_type precision = field_descr->precision(id_field);
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
