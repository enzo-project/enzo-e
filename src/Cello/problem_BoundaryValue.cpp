// See LICENSE_CELLO file for license and copyright information

/// @file     problem_BoundaryValue.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-04-02
/// @brief    Implementation of the default BoundaryValue boundary value class

#include "problem.hpp"

//----------------------------------------------------------------------

void BoundaryValue::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change
  Boundary::pup(p); 
  TRACEPUP; 

  int has_value = (value_!=NULL);
  p | has_value;
  if (has_value){
    if (p.isUnpacking()){
      value_ = new Value;
    }
    p | *value_;
  } else {
    value_ = NULL;
  }
  p | field_list_;
}

//----------------------------------------------------------------------

void BoundaryValue::enforce 
(Block * block, face_enum face, axis_enum axis) const throw()
{
  if ( ! applies_(axis,face)) return;

  if (face == face_all) {
    enforce(block,face_lower,axis);
    enforce(block,face_upper,axis);
  } else if (axis == axis_all) {
    enforce(block,face,axis_x);
    enforce(block,face,axis_y);
    enforce(block,face,axis_z);
  } else {

    Data * data = block->data();
    Field field = data->field();

    if ( ! field.ghosts_allocated() ) {
      ERROR("EnzoBoundary::enforce",
	    "Function called with ghosts not allocated");
    }

    double xm,ym,zm;
    double xp,yp,zp;
    data -> lower(&xm,&ym,&zm);
    data -> upper(&xp,&yp,&zp);

    double t = block->time();

    for (size_t index = 0; index < field_list_.size(); index++) {

      int nx,ny,nz;
      field.size(&nx,&ny,&nz);

      int index_field = field.field_id(field_list_[index]);
      int gx,gy,gz;
      field.ghost_depth(index_field,&gx,&gy,&gz);

      int cx,cy,cz;
      field.centering(index_field, &cx,&cy,&cz);

      int ndx=nx+2*gx+cx;
      int ndy=ny+2*gy+cy;
      int ndz=nz+2*gz+cz;

      double * x = new double [ndx];
      double * y = new double [ndy];
      double * z = new double [ndz];

      data->field_cell_faces(x,y,z,gx,gy,gz,cx,cy,cz);

      void * array = field.values(index_field);

      precision_type precision = field.precision(index_field);

      int ix0=0 ,iy0=0,iz0=0;

      nx = ndx;
      ny = ndy;
      nz = ndz;

      if (axis == axis_x) nx=gx;
      if (axis == axis_y) ny=gy;
      if (axis == axis_z) nz=gz;

      if (face == face_upper) {
	if (axis == axis_x) ix0 = ndx - gx;
	if (axis == axis_y) iy0 = ndy - gy;
	if (axis == axis_z) iz0 = ndz - gz;
      }

      int i0=ix0 + ndx*(iy0 + ndy*iz0);

      bool * mask = 0;

      if (mask_ != nullptr) mask = new bool [nx*ny*nz];

      switch (precision) {
      case precision_single:
	{
	  float * temp = 0;
	  if (mask_ != nullptr) {
	    temp = (float *)array;
	    array = new float [ndx*ndy*ndz];
	  }
	  
	  value_->evaluate((float *)array+i0, t, 
			   ndx,nx,x+ix0, 
			   ndy,ny,y+iy0,
			   ndz,nz,z+iz0);
	  if (mask_ != nullptr) {
	    for (int i=0; i<ndx*ndy*ndz; i++) ((float *)temp)[i]=((float *)array)[i];
	    delete [] ((float*)array);
	    array = temp;
	  }
	}
       	break;
      case precision_double:
	{
	  double * temp = 0;
	  if (mask_ != nullptr) {
	    temp = (double *)array;
	    array = new double [ndx*ndy*ndz];
	  }
	  value_->evaluate((double *)array+i0, t, 
			   ndx,nx,x+ix0, 
			   ndy,ny,y+iy0,
			   ndz,nz,z+iz0);
	  if (mask_ != nullptr) {
	    for (int i=0; i<ndx*ndy*ndz; i++) ((double *)temp)[i]=((double *)array)[i];
	    delete [] ((double *)array);
	    array = temp;
	  }
	}
       	break;
      case precision_extended80:
      case precision_extended96:
      case precision_quadruple:
	{
	  long double * temp = 0;
	  if (mask_ != nullptr) {
	    temp = (long double *)array;
	    array = new long double [ndx*ndy*ndz];
	  }
	  value_->evaluate((long double *)array+i0, t, 
			   ndx,nx,x+ix0, 
			   ndy,ny,y+iy0,
			   ndz,nz,z+iz0);
	  if (mask_ != nullptr) {
	    for (int i=0; i<ndx*ndy*ndz; i++) 
	      ((long double *)temp)[i]=((long double *)array)[i];
	    delete [] ((long double *)array);
	    array = temp;
	  }
	}
       	break;
      }

      delete [] x;
      delete [] y;
      delete [] z;
      delete [] mask;
    }
  }
}

//----------------------------------------------------------------------

template <class T>
void BoundaryValue::copy_(T * field, double * value,
			  int ndx, int ndy, int ndz,
			  int nx,  int ny,  int nz,
			  int ix0, int iy0, int iz0) const throw()
{
  for (int ix=ix0; ix<ix0+nx; ix++) {
    for (int iy=iy0; iy<iy0+ny; iy++) {
      for (int iz=iz0; iz<iz0+nz; iz++) {
	int iv = (ix-ix0) + nx*((iy-iy0) + ny*(iz-iz0));
	int ib = ix + ndx*(iy + ndy*(iz));
	field[ib] = (T) value[iv];
      }
    }
  }
}

//----------------------------------------------------------------------
