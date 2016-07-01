// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeAcceleration.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoComputeAcceleration class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputeAcceleration::EnzoComputeAcceleration 
(FieldDescr * field_descr,
 int         rank,
 int         order)
  : Compute(),
    rank_(rank),
    order_(order)
{
  i_ax_ = (rank_ >= 1) ? field_descr->field_id("acceleration_x") : -1;
  i_ay_ = (rank_ >= 2) ? field_descr->field_id("acceleration_y") : -1;
  i_az_ = (rank_ >= 3) ? field_descr->field_id("acceleration_z") : -1;

  i_p_ =  field_descr->field_id("potential");

  if (order_ != 2 && order_ != 4) {
    ERROR1("EnzoComputeAcceleration",
	   "Unknown order %d", order_);
  }

}

//----------------------------------------------------------------------

void EnzoComputeAcceleration::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Compute::pup(p);

  p | rank_;
  p | order_;
  p | i_ax_;
  p | i_ay_;
  p | i_az_;
  p | i_p_;

}

//----------------------------------------------------------------------

void EnzoComputeAcceleration::compute ( Block * block) throw()
{
  Field field = block->data()->field();
  if (field.precision(0) == precision_single) {
    compute_<float>(block);
  } else if (field.precision(0) == precision_double) {
    compute_<double>(block);
  } else if (field.precision(0) == precision_quadruple) {
    compute_<long double>(block);
  }
}

//----------------------------------------------------------------------

template <typename T>
void EnzoComputeAcceleration::compute_(Block * block)
{
  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  Field field = enzo_block->data()->field();

  T * ax = (rank_ >= 1) ? (T*) field.values(i_ax_) : NULL;
  T * ay = (rank_ >= 2) ? (T*) field.values(i_ay_) : NULL;
  T * az = (rank_ >= 3) ? (T*) field.values(i_az_) : NULL;
  T * p  = (T*) field.values(i_p_);

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);

  int dx,dy,dz;
  dx = 1; 
  dy = mx;
  dz = mx*my;

  int dx2,dy2,dz2;
  dx2 = 2*dx;
  dy2 = 2*dy;
  dz2 = 2*dz;

  double xm,ym,zm;
  double xp,yp,zp;
  double hx,hy,hz;

  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,
		   ym,yp,&hy,
		   zm,zp,&hz);

  if (order_ == 2) {

    if (rank_ == 1) {

      const T fx = 1.0 / (2.0*hx);
      for (int ix=1; ix<mx-1; ix++) {
	int i=ix;
	ax[i] = fx*(p[i+dx] - p[i-dx]);
      }

    } else if (rank_ == 2) {

      const T fx = 1.0 / (2.0*hx);
      const T fy = 1.0 / (2.0*hy);
      for (int iy=1; iy<my-1; iy++) {
	for (int ix=1; ix<mx-1; ix++) {
	  int i=ix + mx*iy;
	  ax[i] = fx*(p[i+dx] - p[i-dx]);
	  ay[i] = fy*(p[i+dy] - p[i-dy]);
	}
      }

    } else if (rank_ == 3) {

      const T fx = 1.0 / (2.0*hx);
      const T fy = 1.0 / (2.0*hy);
      const T fz = 1.0 / (2.0*hz);
      for (int iz=1; iz<mz-1; iz++) {
	for (int iy=1; iy<my-1; iy++) {
	  for (int ix=1; ix<mx-1; ix++) {
	    int i=ix + mx*(iy + my*iz);
	    ax[i] = fx*(p[i+dx] - p[i-dx]);
	    ay[i] = fy*(p[i+dy] - p[i-dy]);
	    az[i] = fz*(p[i+dz] - p[i-dz]);
	  }
	}
      }
    }

  } else if (order_ == 4) {

    if (rank_ == 1) {

      const T fx = 1.0 / (12.0*hx);

      for (int ix=2; ix<mx-2; ix++) {
	int i=ix;
	ax[i] = fx*( -p[i+dx2] + 8*p[i+dx] - 8*p[i-dx] + p[i-dx2]);
      }

    } else if (rank_ == 2) {

      const T fx = 1.0 / (12.0*hx);
      const T fy = 1.0 / (12.0*hy);
      for (int iy=2; iy<my-2; iy++) {
	for (int ix=2; ix<mx-2; ix++) {
	  int i=ix + mx*iy;
	  ax[i] = fx*( -p[i+dx2] + 8*p[i+dx] - 8*p[i-dx] + p[i-dx2]);
	  ay[i] = fy*( -p[i+dy2] + 8*p[i+dy] - 8*p[i-dy] + p[i-dy2]);
	}
      }

    } else if (rank_ == 3) {

      const T fx = 1.0 / (12.0*hx);
      const T fy = 1.0 / (12.0*hy);
      const T fz = 1.0 / (12.0*hz);
      for (int iz=2; iz<mz-2; iz++) {
	for (int iy=2; iy<my-2; iy++) {
	  for (int ix=2; ix<mx-2; ix++) {
	    int i=ix + mx*(iy + my*iz);
	    ax[i] = fx*( -p[i+dx2] + 8*p[i+dx] - 8*p[i-dx] + p[i-dx2]);
	    ay[i] = fy*( -p[i+dy2] + 8*p[i+dy] - 8*p[i-dy] + p[i-dy2]);
	    az[i] = fz*( -p[i+dz2] + 8*p[i+dz] - 8*p[i-dz] + p[i-dz2]);
	  }
	}
      }

    }

  } else {
    ERROR1("EnzoComputeAcceleration",
	   "Unknown order %d", order_);
  }
}

