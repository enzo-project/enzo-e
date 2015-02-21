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
 bool        symmetric,
 int         order)
  : Compute(),
    rank_(rank),
    symmetric_(symmetric),
    order_(order)
{
  i_ax_ = field_descr->field_id("acceleration_x");
  i_ay_ = field_descr->field_id("acceleration_y");
  i_az_ = field_descr->field_id("acceleration_z");

  i_p_ =  field_descr->field_id("potential");

  if (order_ != 2 && order_ != 4 && order_ != 6) {
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
  p | symmetric_;
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

  T * ax = (T*) field.values(i_ax_);
  T * ay = (T*) field.values(i_ay_);
  T * az = (T*) field.values(i_az_);
  T * p  = (T*) field.values(i_p_);

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  int gx,gy,gz;
  field.ghosts (0,&gx,&gy,&gz);

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);

  int dx,dy,dz;
  dx = 1; 
  dy = mx;
  dz = mx*my;

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
      if (symmetric_) {
	const T fx = 1.0 / (2.0*hx);
	for (int ix=1; ix<mx-1; ix++) {
	  int i=ix;
	  ax[i] = fx*(p[i+dx] - p[i-dx]);
	}
      } else { // ! symmetric_
	const T fx = 1.0 / hx;
	for (int ix=0; ix<mx-1; ix++) {
	  int i=ix;
	  ax[i] = fx*(p[i+dx] - p[i]);
	}
      }
    } else if (rank_ == 2) {
      if (symmetric_) {
	const T fx = 1.0 / (2.0*hx);
	const T fy = 1.0 / (2.0*hy);
	for (int iy=1; iy<my-1; iy++) {
	  for (int ix=1; ix<mx-1; ix++) {
	    int i=ix + mx*iy;
	    ax[i] = fx*(p[i+dx] - p[i-dx]);
	    ay[i] = fy*(p[i+dy] - p[i-dy]);
	  }
	}
      } else { // ! symmetric_
	const T fx = 1.0 / hx;
	const T fy = 1.0 / hy;
	for (int iy=0; iy<my-1; iy++) {
	  for (int ix=0; ix<mx-1; ix++) {
	    int i=ix + mx*iy;
	    ax[i] = fx*(p[i+dx] - p[i]);
	    ay[i] = fy*(p[i+dy] - p[i]);
	  }
	}
      }
    } else if (rank_ == 3) {
      if (symmetric_) {
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
      } else { // ! symmetric_
	const T fx = 1.0 / hx;
	const T fy = 1.0 / hy;
	const T fz = 1.0 / hz;
	for (int iz=0; iz<mz-1; iz++) {
	  for (int iy=0; iy<my-1; iy++) {
	    for (int ix=0; ix<mx-1; ix++) {
	      int i=ix + mx*(iy + my*iz);
	      ax[i] = fx*(p[i+dx] - p[i]);
	      ay[i] = fy*(p[i+dy] - p[i]);
	      az[i] = fz*(p[i+dz] - p[i]);
	    }
	  }
	}
      }
    }
  } else {
    ERROR1("EnzoComputeAcceleration",
	   "Unknown order %d", order_);
  }
}

