// See LICENSE_CELLO file for license and copyright information

/// @file     data_Data.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-03-10
/// @brief    Implementation of Data class

#include "data.hpp"

//----------------------------------------------------------------------

Data::Data(int nx, int ny, int nz,
	   int num_field_data,
	   double xm, double xp,
	   double ym, double yp,
	   double zm, double zp,
	   FieldDescr * field_descr,
	   ParticleDescr * particle_descr) throw ()
  : num_field_data_(num_field_data),
    field_data_(),
    particle_data_(),
    flux_data_()
{
  if (field_descr == nullptr)
    field_descr = cello::field_descr();
  if (particle_descr == nullptr)
    particle_descr = cello::particle_descr();
  
  // Initialize field data
  field_data_.resize(num_field_data);
  for (size_t i=0; i<field_data_.size(); i++) {
    field_data_[i] = new FieldData (field_descr,nx,ny,nz);
  }
  // Initialize particle data
  particle_data_ = new ParticleData;
  particle_data_->allocate(particle_descr);

  // Initialize flux data
  flux_data_ = new FluxData;
  
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
    field_data_[i] = nullptr;
  }
  num_field_data_ = 0;

  delete particle_data_;
  particle_data_ = nullptr;

  delete flux_data_;
  flux_data_ = nullptr;
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

void Data::field_cell_faces (double * x, double * y, double * z,
			     int gx, int gy, int gz,
			     int cx, int cy, int cz) const
{
  // cx, cy, and cz are expected to be 1 or 0. If its 0, then the
  // position along the corresponding dimension are for cell centers
  double hx,hy,hz;
  field_cell_width (&hx,&hy,&hz);

  double xm,ym,zm;
  this->lower(&xm,&ym,&zm);

  int nx,ny,nz;
  field_data()->size(&nx,&ny,&nz);

  int ixm = -gx;
  int iym = -gy;
  int izm = -gz;

  int ixp = nx+gx+cx;
  int iyp = ny+gy+cy;
  int izp = nz+gz+cz;

  double dx = (cx == 0) ? 0.5 : 0;
  double dy = (cy == 0) ? 0.5 : 0;
  double dz = (cz == 0) ? 0.5 : 0;
  
  for (int ix=ixm; ix<ixp; ix++) x[ix-ixm] = xm + (ix+dx)*hx;
  for (int iy=iym; iy<iyp; iy++) y[iy-iym] = ym + (iy+dy)*hy;
  for (int iz=izm; iz<izp; iz++) z[iz-izm] = zm + (iz+dz)*hz;
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

void Data::allocate () throw()
{
  // allocate Block scalar storage
  scalar_data_long_double_.allocate(cello::scalar_descr_long_double());
  scalar_data_double_     .allocate(cello::scalar_descr_double());
  scalar_data_int_        .allocate(cello::scalar_descr_int());
  scalar_data_sync_       .allocate(cello::scalar_descr_sync());
  scalar_data_void_       .allocate(cello::scalar_descr_void());
  // allocate Block Field storage
  for (size_t i=0; i<field_data_.size(); i++) {
    field_data_[i]->set_history_(cello::field_descr());
    field_data_[i]->allocate_permanent(cello::field_descr(),true);
  }
  // initialize Flux field list
}

//======================================================================

void Data::copy_(const Data & data) throw()
{
  num_field_data_ = data.num_field_data_;
  field_data_.resize(data.field_data_.size());
  for (size_t i=0; i<field_data_.size(); i++) {
    field_data_[i] = new FieldData (*(data.field_data_[i]));
  }
  particle_data_ = new ParticleData (*data.particle_data_);
  flux_data_ = new FluxData (*data.flux_data_);
}
