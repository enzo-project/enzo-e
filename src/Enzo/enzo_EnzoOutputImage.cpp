// $Id: enzo_EnzoOutputImage.cpp 2093 2011-03-12 01:17:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoOutputImage.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
/// @brief    Implementation of the EnzoOutputImage class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoOutputImage::EnzoOutputImage() throw ()
  : Output()
{
}

//----------------------------------------------------------------------

EnzoOutputImage::~EnzoOutputImage() throw ()
{
}

//======================================================================

void EnzoOutputImage::write
(
 int index,
 Mesh * mesh,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
  ) const throw()
{
  if (mesh->dimension() != 2) {
    WARNING("EnzoOutputImage::write",
	    "EnzoOutputImage only supports 2D problems");
  }

  // Get field_descr (eventually)

  Block * block = mesh->root_patch()->local_block(0);
  if (block == NULL) return;
  FieldBlock       * field_block = block->field_block();
  const FieldDescr * field_descr = field_block->field_descr();

  Monitor * monitor = Monitor::instance();

  // Open file if writing a single block
  if (root_call) {

    // Get file name
    std::string file_prefix = expand_file_name (cycle,time);
    std::string field_name  = field_descr->field_name(index);
    std::string file_name = file_prefix + "-" + field_name + ".png";

    // Get mesh size
    int nxm,nym,nzm;
    mesh->root_patch()->size (&nxm, &nym, &nzm);
    // Create image 
    monitor->image_open(file_name,nxm,nym);
  }

  ItPatch it_patch (mesh);
  while (Patch * patch = ++it_patch) {

    // Write patch contribution 
    // NO OFFSET: ASSUMES ROOT PATCH
    write (index, patch, mesh, cycle,time,false,
	   0,0,0);
  }

  if (root_call) {
    // Complete geterating image if writing a single block
    double min=0.0;
    double max=1.0;
    if (root_call) monitor->image_close(min,max);
  }

}

//----------------------------------------------------------------------

void EnzoOutputImage::write
(
 int index,
 Patch * patch,
 Mesh  * mesh,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
 ) const throw()
{
  if (mesh->dimension() != 2) {
    WARNING("EnzoOutputImage::write",
	    "EnzoOutputImage only supports 2D problems");
  }

  Block * block = patch->local_block(0);
  FieldBlock       * field_block = block->field_block();
  const FieldDescr * field_descr = field_block->field_descr();

  Monitor * monitor = Monitor::instance();

  // Open file if writing a single block
  if (root_call) {

    // Get file name
    std::string file_prefix = expand_file_name (cycle,time);
    std::string field_name = field_descr->field_name(index);
    std::string file_name = file_prefix + "-" + field_name + ".png";

    // Get patch size
    int nxp, nyp;
    patch->size (&nxp, &nyp);

    // Create image 
    monitor->image_open(file_name,nxp,nyp);
  }

  ItBlock it_block (patch);
  while (Block * block = ++it_block) {
    // Get block size
    int nxb,nyb,nzb;
    FieldBlock * field_block = block->field_block();
    field_block->size(&nxb,&nyb,&nzb);

    int ix,iy,iz;
    block->index_patch(&ix,&iy,&iz);

    write (index, block, patch, mesh, cycle,time,false,
	   ix0+ix*nxb,
	   iy0+iy*nyb,
	   iz0+iz*nzb);
  }

  if (root_call) {
    // Complete geterating image if writing a single block
    double min=0.0;
    double max=1.0;
    if (root_call) monitor->image_close(min,max);
  }

}

//----------------------------------------------------------------------

void EnzoOutputImage::write
(
 int index,
 Block * block,
 Patch * patch,
 Mesh  * mesh,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
) const throw()
{
  FieldBlock       * field_block = block->field_block();
  const FieldDescr * field_descr = field_block->field_descr();

  // Get block size
  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  std::string file_prefix = expand_file_name (cycle,time);

  Monitor * monitor = Monitor::instance();

  // Get ghost depth
  int gx,gy,gz;
  field_descr->ghosts(index,&gx,&gy,&gz);

  // Get array dimensions
  int ndx,ndy,ndz;
  ndx=nx+2*gx;
  ndy=ny+2*gy;
  ndz=nz+2*gz;

  if (root_call) {
    // Open file if writing a single block
    std::string field_name = field_descr->field_name(index);
    std::string file_name = file_prefix + "-" + field_name + ".png";
    monitor->image_open(file_name,nx,ny);
  }

  // Add block contribution to image

  Scalar * field_unknowns = (Scalar *)field_block->field_unknowns(index);
    
  monitor->image_reduce (field_unknowns, 
			 ndx,ndy,ndz, 
			 nx,ny,nz,
			 ix0,iy0,iz0,
			 axis_z, reduce_sum);

  if (root_call) {
    // Complete geterating image if writing a single block
    double min=0.0;
    double max=1.0;
    if (root_call) monitor->image_close(min,max);
  }
}
