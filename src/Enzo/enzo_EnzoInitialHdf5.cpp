// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitiaHdf5.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-06-23
/// @brief    Read initial conditions from HDF5
///           (multi-scale cosmological initial conditions)

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialHdf5::EnzoInitialHdf5
(int cycle,
 double time,
 const EnzoConfig * enzo_config) throw()
  : Initial(cycle,time),
    field_files_     (enzo_config->initial_hdf5_field_files),
    field_datasets_  (enzo_config->initial_hdf5_field_datasets),
    field_coords_    (enzo_config->initial_hdf5_field_coords),
    field_names_     (enzo_config->initial_hdf5_field_names),
    particle_files_     (enzo_config->initial_hdf5_particle_files),
    particle_datasets_  (enzo_config->initial_hdf5_particle_datasets),
    particle_types_     (enzo_config->initial_hdf5_particle_types),
    particle_attributes_(enzo_config->initial_hdf5_particle_attributes)
{
}

//----------------------------------------------------------------------

EnzoInitialHdf5::EnzoInitialHdf5() throw ()
{
  INCOMPLETE("EnzoInitialHdf5::EnzoInitialHdf5");
}

//----------------------------------------------------------------------

EnzoInitialHdf5::~EnzoInitialHdf5() throw ()
{
  INCOMPLETE("EnzoInitialHdf5::~EnzoInitialHdf5");
}

//----------------------------------------------------------------------

void EnzoInitialHdf5::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
  p | field_files_;
  p | field_datasets_;
  p | field_coords_;
  p | field_names_;
  p | particle_files_;
  p | particle_datasets_;
  p | particle_types_;
  p | particle_attributes_;
}

void EnzoInitialHdf5::enforce_block
  ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
    ) throw()
{

  // Get the domain extents

  double lower_domain[3];
  double upper_domain[3];
  hierarchy->lower(lower_domain,
		   lower_domain+1,
		   lower_domain+2);
  hierarchy->upper(upper_domain,
		   upper_domain+1,
		   upper_domain+2);

  double lower_block[3];
  double upper_block[3];

  block->lower(lower_block,
	       lower_block+1,
	       lower_block+2);

  block->upper(upper_block,
	       upper_block+1,
	       upper_block+2);

  for (int i=0; i<field_files_.size(); i++) {

    std::string file_name = field_files_[i];

    Field field = block->data()->field();
    enzo_float * data = (enzo_float *) field.values(field_names_[i]);

    // Block size

    int mx,my,mz;
    int nx,ny,nz;
    int gx,gy,gz;
 
    field.dimensions (0,&mx,&my,&mz);
    field.size         (&nx,&ny,&nz);
    field.ghost_depth(0,&gx,&gy,&gz);

    // Open the field file

    FileHdf5 file("./",file_name);

    file.file_open();

    // create memory space
    file.mem_create (mx,my,mz,nx,ny,nz,gx,gy,gz);

    // Read the domain dimensions

    const int ix = field_coords_[i].find ("x");
    const int iy = field_coords_[i].find ("y");
    const int iz = field_coords_[i].find ("z");

    ASSERT3 ("EnzoInitialHdf5::enforce_block()",
	    "bad field coordinates %d %d %d",
	     ix,iy,iz,
	     ((ix<4)&&(iy<4)&&(iz<4)) &&
	     ((ix != iy) || (iy==-1 && iz == -1)) &&
	     ((ix != iy && iy != iz) || (iz == -1)));
    
    int m4[4] = {0};

    int type = type_unknown;
    file. data_open (field_datasets_[i], &type,
		     m4,m4+1,m4+2,m4+3);

    mx = m4[ix];
    my = m4[iy];
    mz = m4[iz];
    
    // compute cell widths
    double h4[4] = {1};
    h4[ix] = (upper_domain[0] - lower_domain[0]) / m4[ix];
    h4[iy] = (upper_domain[1] - lower_domain[1]) / m4[iy];
    h4[iz] = (upper_domain[2] - lower_domain[2]) / m4[iz];

    // determine offsets
    
    int o4[4] = {0};
    o4[ix] = (lower_block[0] - lower_domain[0]) / h4[ix];
    o4[iy] = (lower_block[1] - lower_domain[1]) / h4[iy];
    o4[iz] = (lower_block[2] - lower_domain[2]) / h4[iz];

    int n4[4] = {1};
    n4[ix] = (upper_block[0] - lower_block[0]) / h4[ix];
    n4[iy] = (upper_block[1] - lower_block[1]) / h4[iy];
    n4[iz] = (upper_block[2] - lower_block[2]) / h4[iz];
      
    // domain |------------------|
  // block       |-----|
  // offset |....|
  // length      |.....|
  //
  // assert block length == expected block length (for now)
  // hx = (upper_domain - lower_domain) / m4
  // offset = (lower_block - lower_domain)/hx
  // length = (upper_block - lower_block) / hx

  

    // open the dataspace
    file. data_slice
      (m4[0],m4[1],m4[2],m4[3],
       n4[0],n4[1],n4[2],n4[3],
       o4[0],o4[1],o4[2],o4[3]);

    // input domain size
  
    file.data_read (data);
    file.data_close();
    file.file_close();

  }

  for (int i=0; i<particle_files_.size(); i++) {
    CkPrintf ("DEBUG_INITIAL particle file %s\n",
	      particle_files_[i].c_str());
  }  
  fflush(stdout);

}

//======================================================================

