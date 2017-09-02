// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitiaHdf5.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-06-23
/// @brief    Read initial conditions from HDF5
///           (multi-scale cosmological initial conditions)

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialMusic::EnzoInitialMusic
(int cycle,
 double time,
 const EnzoConfig * enzo_config) throw()
  : Initial(cycle,time),
    level_(enzo_config->mesh_max_level),
    field_files_     (enzo_config->initial_music_field_files),
    field_datasets_  (enzo_config->initial_music_field_datasets),
    field_coords_    (enzo_config->initial_music_field_coords),
    field_names_     (enzo_config->initial_music_field_names),
    particle_files_     (enzo_config->initial_music_particle_files),
    particle_datasets_  (enzo_config->initial_music_particle_datasets),
    particle_coords_    (enzo_config->initial_music_particle_coords),
    particle_types_     (enzo_config->initial_music_particle_types),
    particle_attributes_(enzo_config->initial_music_particle_attributes)
{ }

//----------------------------------------------------------------------

void EnzoInitialMusic::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change

  p | level_;
  
  p | field_files_;
  p | field_datasets_;
  p | field_coords_;
  p | field_names_;

  p | particle_files_;
  p | particle_datasets_;
  p | particle_coords_;
  p | particle_types_;
  p | particle_attributes_;
}

void EnzoInitialMusic::enforce_block
  ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
    ) throw()
{

  if (block->level() != level_) return;
  // Get the domain extents

  double lower_domain[3];
  double upper_domain[3];
  hierarchy->lower(lower_domain, lower_domain+1, lower_domain+2);
  hierarchy->upper(upper_domain, upper_domain+1, upper_domain+2);

  double lower_block[3];
  double upper_block[3];
  block->lower(lower_block, lower_block+1, lower_block+2);
  block->upper(upper_block, upper_block+1, upper_block+2);

  Field field = block->data()->field();

  for (int i=0; i<field_files_.size(); i++) {

    std::string file_name = field_files_[i];

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

    // Read the domain dimensions

    const int IX = field_coords_[i].find ("x");
    const int IY = field_coords_[i].find ("y");
    const int IZ = field_coords_[i].find ("z");

    ASSERT3 ("EnzoInitialMusic::enforce_block()",
	    "bad field coordinates %d %d %d",
	     IX,IY,IZ,
	     ((IX<4)&&(IY<4)&&(IZ<4)) &&
	     ((IX != IY) || (IY==-1 && IZ == -1)) &&
	     ((IX != IY && IY != IZ) || (IZ == -1)));
    
    int m4[4] = {0};
    int type = type_unknown;
    file. data_open (field_datasets_[i], &type,
		     m4,m4+1,m4+2,m4+3);

    // compute cell widths
    double h4[4] = {1};
    h4[IX] = (upper_domain[0] - lower_domain[0]) / m4[IX];
    h4[IY] = (upper_domain[1] - lower_domain[1]) / m4[IY];
    h4[IZ] = (upper_domain[2] - lower_domain[2]) / m4[IZ];

    // determine offsets
    int o4[4] = {0};
    o4[IX] = (lower_block[0] - lower_domain[0]) / h4[IX];
    o4[IY] = (lower_block[1] - lower_domain[1]) / h4[IY];
    o4[IZ] = (lower_block[2] - lower_domain[2]) / h4[IZ];

    int n4[4] = {1};
    n4[IX] = (upper_block[0] - lower_block[0]) / h4[IX];
    n4[IY] = (upper_block[1] - lower_block[1]) / h4[IY];
    n4[IZ] = (upper_block[2] - lower_block[2]) / h4[IZ];
      
    // open the dataspace
    file. data_slice
      (m4[0],m4[1],m4[2],m4[3],
       n4[0],n4[1],n4[2],n4[3],
       o4[0],o4[1],o4[2],o4[3]);

    // create memory space
    file.mem_create (n4[IX],n4[IY],n4[IZ],
		     n4[IX],n4[IY],n4[IZ],
		     0,0,0);
    
    // input domain size
  
    enzo_float * data = new enzo_float[nx*ny*nz];

    file.data_read (data);

    enzo_float * array = (enzo_float *) field.values(field_names_[i]);

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  array[ix+gx+mx*(iy+gy+my*(iz+gz))] = data[ix+n4[IX]*(iy+n4[IY]*iz)];
	}
      }
    }

    delete [] data;
    
    file.data_close();
    file.file_close();

  }


  for (int i=0; i<particle_files_.size(); i++) {

    std::string file_name = particle_files_[i];

    // Open the particle file

    FileHdf5 file("./",file_name);

    file.file_open();

    // Open the dataset
    int m4[4] = {0};
    int type = type_unknown;
    file. data_open (particle_datasets_[i], &type,
		     m4,m4+1,m4+2,m4+3);

    // coordinate mapping
    const int IX = particle_coords_[i].find ("x");
    const int IY = particle_coords_[i].find ("y");
    const int IZ = particle_coords_[i].find ("z");


    // compute cell widths
    double h4[4] = {1};
    h4[IX] = (upper_domain[0] - lower_domain[0]) / m4[IX];
    h4[IY] = (upper_domain[1] - lower_domain[1]) / m4[IY];
    h4[IZ] = (upper_domain[2] - lower_domain[2]) / m4[IZ];

    // determine offsets
    int o4[4] = {0};
    o4[IX] = (lower_block[0] - lower_domain[0]) / h4[IX];
    o4[IY] = (lower_block[1] - lower_domain[1]) / h4[IY];
    o4[IZ] = (lower_block[2] - lower_domain[2]) / h4[IZ];

    int n4[4] = {1};
    n4[IX] = (upper_block[0] - lower_block[0]) / h4[IX];
    n4[IY] = (upper_block[1] - lower_block[1]) / h4[IY];
    n4[IZ] = (upper_block[2] - lower_block[2]) / h4[IZ];

    // open the dataspace
    file. data_slice
      (m4[0],m4[1],m4[2],m4[3],
       n4[0],n4[1],n4[2],n4[3],
       o4[0],o4[1],o4[2],o4[3]);

    // create memory space

    int nx,ny,nz;
    field.size (&nx,&ny,&nz);
 
    file.mem_create (nx,ny,nz,nx,ny,nz,0,0,0);

    enzo_float * data = new enzo_float[nx*ny*nz];

    // read data and close file
    file.data_read (data);
    file.data_close();
    file.file_close();

    // Create particles and initialize them

    Particle particle = block->data()->particle();

    const int it = particle.type_index(particle_types_[i]);
    const int ia = particle.attribute_index(it,particle_attributes_[i]);
    const int da = particle.stride(it,ia);

    const int np = nx*ny*nz;

    // insert particles if they don't exist yet
    if (particle.num_particles(it) == 0) {
      particle.insert_particles(it,np);
      block->simulation()->monitor_insert_particles(np);
    }

    // read particle attribute
    for (int ip=0; ip<np; ip++) {
      int ib,io;
      particle.index(ip,&ib,&io);
      enzo_float * array = (enzo_float *) particle.attribute_array(it,ia,ib);
      array[io] = data[ip];
    }

    // update positions with displacements
    if (particle_datasets_[i] == "ParticleDisplacements_x") {
      for (int iz=0; iz<nz; iz++) {
	for (int iy=0; iy<ny; iy++) {
	  for (int ix=0; ix<nx; ix++) {
	    int ip = ix + nx*(iy + ny*iz);
	    int ib,io;
	    particle.index(ip,&ib,&io);
	    enzo_float * array = (enzo_float *)
	      particle.attribute_array(it,ia,ib);
	    array[io] += lower_block[0] + (ix+0.5)*h4[IX];
	  }
	}
      }
    } else if (particle_datasets_[i] == "ParticleDisplacements_y") {
      for (int iz=0; iz<nz; iz++) {
	for (int iy=0; iy<ny; iy++) {
	  for (int ix=0; ix<nx; ix++) {
	    int ip = ix + nx*(iy + ny*iz);
	    int ib,io;
	    particle.index(ip,&ib,&io);
	    enzo_float * array = (enzo_float *)
	      particle.attribute_array(it,ia,ib);
	    array[io] += lower_block[1] + (iy+0.5)*h4[IY];
	  }
	}
      }
    } else if (particle_datasets_[i] == "ParticleDisplacements_z") {
      for (int iz=0; iz<nz; iz++) {
	for (int iy=0; iy<ny; iy++) {
	  for (int ix=0; ix<nx; ix++) {
	    int ip = ix + nx*(iy + ny*iz);
	    int ib,io;
	    particle.index(ip,&ib,&io);
	    enzo_float * array = (enzo_float *)
	      particle.attribute_array(it,ia,ib);
	    array[io] += lower_block[2] + (iz+0.5)*h4[IZ];
	  }
	}
      }
    }
    
    delete [] data;
    data = NULL;
  }  
}
