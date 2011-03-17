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
  INCOMPLETE("EnzoOutputImage::EnzoOutputImage");
}

//----------------------------------------------------------------------

EnzoOutputImage::~EnzoOutputImage() throw ()
{
  INCOMPLETE("EnzoOutputImage::!EnzoOutputImage");
}

//======================================================================

void EnzoOutputImage::write
(
 Mesh * mesh,
 int cycle,
 double time,
 bool root_call
  ) const throw()
{
  // Check that dimension == 2

  INCOMPLETE("EnzoOutputImage::write");

  if (mesh->dimension() != 2) {
    WARNING("EnzoOutputImage::write",
	    "EnzoOutputImage only supports 2D problems");
  }

  if (root_call) {
    
    // create mesh image

    // open file
    
    // write header
  }

  ItPatch it_patch (mesh);
  while (Patch * patch = ++it_patch) {
    write (mesh,patch,cycle,time,false);
  }

  if (root_call) {
    // write footer
    // close file
  }

}

//----------------------------------------------------------------------

void EnzoOutputImage::write
(
 Mesh  * mesh,
 Patch * patch,
 int cycle,
 double time,
 bool root_call
 ) const throw()
{
  INCOMPLETE("EnzoOutputImage::write");

  if (mesh->dimension() != 2) {
    WARNING("EnzoOutputImage::write",
	    "EnzoOutputImage only supports 2D problems");
  }

  if (root_call) {
    
    // create mesh image

    // open file
    
    // write header
  }

  ItBlock it_block (patch);
  while (Block * block = ++it_block) {
    write (mesh, patch, block,cycle,time,false);
  }

  if (root_call) {
    // write footer
    // close file
  }

}

//----------------------------------------------------------------------

void EnzoOutputImage::write
(
 Mesh  * mesh,
 Patch * patch,
 Block * block,
 int cycle,
 double time,
 bool root_call
) const throw()
{
  INCOMPLETE("EnzoOutputImage::write");

  printf ("write(block) file_name = %s\n",
	  expand_file_name (cycle,time).c_str());


}

//----------------------------------------------------------------------

// void EnzoSimulationMpi::output_images_
// (
//  Block * block,
//  const char * file_format,
//  int cycle,
//  int cycle_skip
//  ) throw ()
// {

  // if (! (cycle_skip && cycle % cycle_skip == 0)) return;

  // Monitor * monitor = Monitor::instance();
  // FieldBlock *       field_block = block->field_block();
  // const FieldDescr * field_descr = field_block->field_descr();
  // int nx,ny,nz;
  // int gx,gy,gz;
  // int mx,my,mz;
  // field_block->size(&nx,&ny,&nz);
  // int count = field_descr->field_count();
  // for (int index = 0; index < count; index++) {
  //   EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);
  //   field_descr->ghosts(index,&gx,&gy,&gz);
  //   mx=nx+2*gx;
  //   my=ny+2*gy;
  //   mz=nz+2*gz;
  //   char filename[255];
  //   std::string field_name = field_descr->field_name(index);
  //   Scalar * field_values = (Scalar *)field_block->field_values(index);
  //   sprintf (filename,file_format,
  // 	     enzo_block->CycleNumber,index);
  //   monitor->image (filename, field_values, mx,my,mz, 2, reduce_sum, 
  // 		    0.0, 1.0);
  // }
