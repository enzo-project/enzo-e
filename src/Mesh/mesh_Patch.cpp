// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the Patch class

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Patch::Patch
(
 Factory * factory, 
 GroupProcess * group_process,
 int nx,  int ny,  int nz,
 int nbx, int nby, int nbz,
 double xm, double ym, double zm,
 double xp, double yp, double zp
) throw()
  : factory_(factory),
    group_process_(group_process),
    layout_(new Layout (nbx,nby,nbz))
{
  // Check 
  if ( ! ((nx >= nbx) && (ny >= nby) && (nz >= nbz))) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,
	     "Patch size (%d,%d,%d) must be larger than blocking (%d,%d,%d)",
	     nx,ny,nz,nbx,nby,nbz);
    ERROR("Patch::Patch", buffer);
  }

  size_[0] = nx;
  size_[1] = ny;
  size_[2] = nz;

  blocking_[0] = nbx;
  blocking_[1] = nby;
  blocking_[2] = nbz;

  lower_[0] = xm;
  lower_[1] = ym;
  lower_[2] = zm;

  upper_[0] = xp;
  upper_[1] = yp;
  upper_[2] = zp;

}

//----------------------------------------------------------------------

Patch::~Patch() throw()
{
  delete layout_;
}

//----------------------------------------------------------------------

// Patch::Patch(const Patch & patch,
// 	     FieldDescr *  field_descr) throw()
// {
//   deallocate_blocks();

//   allocate_blocks(field_descr);
// }

//----------------------------------------------------------------------

// void Patch::set_lower(double xm, double ym, double zm) throw ()
// {
//   lower_[0] = xm;
//   lower_[1] = ym;
//   lower_[2] = zm;
// };

// //----------------------------------------------------------------------
// void Patch::set_upper(double xp, double yp, double zp) throw ()
// {
//   upper_[0] = xp;
//   upper_[1] = yp;
//   upper_[2] = zp;
// };

// //----------------------------------------------------------------------

// Patch & Patch::operator= (const Patch & patch) throw()
// {
//   deallocate_blocks();
//   allocate_blocks();
//   return *this;
// }

//----------------------------------------------------------------------

void Patch::size (int * npx, int * npy, int * npz) const throw()
{
  if (npx) *npx = size_[0];
  if (npy) *npy = size_[1];
  if (npz) *npz = size_[2];
}

//----------------------------------------------------------------------

void Patch::blocking (int * nbx, int * nby, int * nbz) const throw()
{
  *nbx = blocking_[0];
  *nby = blocking_[1];
  *nbz = blocking_[2];
}

//----------------------------------------------------------------------

Layout * Patch::layout () const throw()
{
  return layout_;
}

//----------------------------------------------------------------------
  
void Patch::lower(double * xm, double * ym, double * zm) const throw ()
{
  if (xm) *xm = lower_[0];
  if (ym) *ym = lower_[1];
  if (zm) *zm = lower_[2];
}

//----------------------------------------------------------------------
void Patch::upper(double * xp, double * yp, double * zp) const throw ()
{
  if (xp) *xp = upper_[0];
  if (yp) *yp = upper_[1];
  if (zp) *zp = upper_[2];
}

//----------------------------------------------------------------------

size_t Patch::num_local_blocks() const  throw()
{
  int rank = group_process_->rank();
  return layout_->local_count(rank);
}

//======================================================================
