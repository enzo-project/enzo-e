// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Data.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-03-10
/// @brief    Implementation of Data class

#include "mesh.hpp"

//----------------------------------------------------------------------

Data::Data(FieldDescr * field_descr,
	     int nx, int ny, int nz,
	     int num_field_data,
	     double xm, double xp,
	     double ym, double yp,
	     double zm, double zp) throw ()
  : num_field_data_(num_field_data),
    field_data_()
{
  // Initialize field_data_[]
  field_data_.resize(num_field_data);
  for (size_t i=0; i<field_data_.size(); i++) {
    field_data_[i] = new FieldData (field_descr,nx,ny,nz);
  }
  lower_[0] = xm;
  lower_[1] = ym;
  lower_[2] = zm;
  upper_[0] = xp;
  upper_[1] = yp;
  upper_[2] = zp;
}

//----------------------------------------------------------------------

Data::~Data() throw ()
{
  // Deallocate field_data_[]
  for (size_t i=0; i<field_data_.size(); i++) {
    delete field_data_[i];
    field_data_[i] = 0;
  }
  num_field_data_ = 0;
}

//----------------------------------------------------------------------

Data::Data(const Data & data) throw ()
/// @param     data  Object being copied
:
  num_field_data_(data.num_field_data_)

{
  copy_(data);
}

//----------------------------------------------------------------------

Data & Data::operator= (const Data & data) throw ()
/// @param     data  Source object of the assignment
/// @return    The target assigned object
{
  copy_(data);
  return *this;
}

//----------------------------------------------------------------------

void Data::allocate (const FieldDescr * field_descr) throw()
{
  for (size_t i=0; i<field_data_.size(); i++) {
    field_data_[i]->allocate_permanent(true);
  }
}

//----------------------------------------------------------------------

void Data::field_cell_width 
(double * hx, 
 double * hy,
 double * hz) const
{
  double xm,ym,zm;
  this->lower(&xm,&ym,&zm);
  double xp,yp,zp;
  this->upper(&xp,&yp,&zp);
  field_data()->cell_width(xm,xp,hx, ym,yp,hy, zm,zp,hz);
}   

//----------------------------------------------------------------------

void Data::field_cells (double * x, double * y, double * z,
			 int gx, int gy, int gz) const
{
  double hx,hy,hz;

  field_cell_width (&hx,&hy,&hz);

  double xm,ym,zm;
  this->lower(&xm,&ym,&zm);

  int nx,ny,nz;
  field_data()->size(&nx,&ny,&nz);

  int ixm = -gx;
  int iym = -gy;
  int izm = -gz;

  int ixp = nx+gx;
  int iyp = ny+gy;
  int izp = nz+gz;
  
  for (int ix=ixm; ix<ixp; ix++) x[ix-ixm] = xm + (ix+0.5)*hx;
  for (int iy=iym; iy<iyp; iy++) y[iy-iym] = ym + (iy+0.5)*hy;
  for (int iz=izm; iz<izp; iz++) z[iz-izm] = zm + (iz+0.5)*hz;
}

//======================================================================

void Data::copy_(const Data & data) throw()
{
  num_field_data_ = data.num_field_data_;
  field_data_.resize(data.field_data_.size());
  for (size_t i=0; i<field_data_.size(); i++) {
    field_data_[i] = new FieldData (*(data.field_data_[i]));
  }
}
