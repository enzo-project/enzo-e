// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitiaHdf5.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-05-13
/// @brief    Read initial conditions from HDF5
#include "enzo.hpp"
#include <chrono>
#include <thread>

//----------------------------------------------------------------------

EnzoInitialHdf5::EnzoInitialHdf5
(int cycle,
 double time,
 const EnzoConfig * enzo_config,
 int level) throw()
  : Initial(cycle,time),
    level_(level),
    field_files_     (enzo_config->initial_hdf5_field_files),
    field_datasets_  (enzo_config->initial_hdf5_field_datasets),
    field_coords_    (enzo_config->initial_hdf5_field_coords),
    field_names_     (enzo_config->initial_hdf5_field_names),
    particle_files_     (enzo_config->initial_hdf5_particle_files),
    particle_datasets_  (enzo_config->initial_hdf5_particle_datasets),
    particle_coords_    (enzo_config->initial_hdf5_particle_coords),
    particle_types_     (enzo_config->initial_hdf5_particle_types),
    particle_attributes_(enzo_config->initial_hdf5_particle_attributes)
{
}

//----------------------------------------------------------------------

void EnzoInitialHdf5::pup (PUP::er &p)
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

//----------------------------------------------------------------------

void EnzoInitialHdf5::enforce_block
( Block * block, const Hierarchy * hierarchy_unused ) throw()
{

  if (block->level() != level_) return;

  // Get the grid size at level_
  double lower_domain[3];
  double upper_domain[3];

  Hierarchy * hierarchy = cello::simulation()->hierarchy();
  
  hierarchy->lower(lower_domain, lower_domain+1, lower_domain+2);
  hierarchy->upper(upper_domain, upper_domain+1, upper_domain+2);

  double lower_block[3];
  double upper_block[3];
  block->lower(lower_block, lower_block+1, lower_block+2);
  block->upper(upper_block, upper_block+1, upper_block+2);

  Field field = block->data()->field();

  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;
  int n4[4],IX,IY,IZ;
  double h4[4];
  
  static std::map<std::string,int> close_count;

  union { void * data; float * data_float; double * data_double; };

  for (size_t index=0; index<field_files_.size(); index++) {

    const int type_data = read_dataset_
      (&data,index,block,
       field_files_[index],
       field_coords_[index],
       field_datasets_[index],
       &mx,&my,&mz,&nx,&ny,&nz,&gx,&gy,&gz,n4,h4,&IX,&IY,&IZ);
    
    enzo_float * array = (enzo_float *) field.values(field_names_[index]);

    if (type_data == type_single) {

      copy_field_data_to_array_
	(array,data_float,mx,my,mz,nx,ny,nz,gx,gy,gz,n4,IX,IY);
      
    } else if (type_data == type_double) {

      copy_field_data_to_array_
	(array,data_double,mx,my,mz,nx,ny,nz,gx,gy,gz,n4,IX,IY);
    }

    if (type_data == type_single) {
      delete [] data_float;
    } else if (type_data == type_double) {
      delete [] data_double;
    }
  }

  for (size_t index=0; index<particle_files_.size(); index++) {

    const int type_data = read_dataset_
      (&data, index, block,
       particle_files_[index],
       particle_coords_[index],
       particle_datasets_[index],
       &mx,&my,&mz,&nx,&ny,&nz,&gx,&gy,&gz,n4,h4,&IX,&IY,&IZ);

    // Create particles and initialize them

    Particle particle = block->data()->particle();

    const int it = particle.type_index(particle_types_[index]);
    const int ia = particle.attribute_index(it,particle_attributes_[index]);

    const int np = nx*ny*nz;

    // insert particles if they don't exist yet
    if (particle.num_particles(it) == 0) {
      particle.insert_particles(it,np);
      enzo::simulation()->data_insert_particles(np);
    }

    // read particle attribute
    union {
      void *   array;
      float *  array_float;
      double * array_double;
    };

    const int type_array = particle.attribute_type(it,ia);

    if (type_array == type_single) {
      if (type_data == type_single) {
	copy_particle_data_to_array_
	  (array_float,data_float,particle,it,ia,np);
      } else if (type_data == type_double) {
	copy_particle_data_to_array_
	  (array_float,data_double,particle,it,ia,np);
      }
    } else if (type_array == type_double) {
      if (type_data == type_single) {
	copy_particle_data_to_array_
	  (array_double,data_float,particle,it,ia,np);
      } else if (type_data == type_double) {
	copy_particle_data_to_array_
	  (array_double,data_double,particle,it,ia,np);
      }
    } else {
      ERROR3 ("EnzoInitialHdf5::enforce_block()",
	      "Unsupported particle precision %s for "
	      "particle type %s attribute %s",
	      cello::precision_name[type_array],
	      particle.type_name(it).c_str(),
	      particle.attribute_name(it,ia).c_str());
    }
    
    if (type_data == type_single) {
      delete [] data_float;
    } else if (type_data == type_double) {
      delete [] data_double;
    }

    data = NULL;

    // update positions with displacements
    if (type_array == type_single) {

      if (particle_datasets_[index] == "ParticleDisplacements_x") 
        update_particle_displacements_
          (array_float,nx,ny,nz,
           particle,it,ia,lower_block[0],h4[IX],0);
      else if (particle_datasets_[index] == "ParticleDisplacements_y") 
        update_particle_displacements_
          (array_float,nx,ny,nz,
           particle,it,ia,lower_block[1],h4[IY],1);
      else if (particle_datasets_[index] == "ParticleDisplacements_z") 
        update_particle_displacements_
          (array_float,nx,ny,nz,
           particle,it,ia,lower_block[2],h4[IZ],2);

    } else { // (type_array != type_single) {

      if (particle_datasets_[index] == "ParticleDisplacements_x") 
        update_particle_displacements_
          (array_double,nx,ny,nz,
           particle,it,ia,lower_block[0],h4[IX],0);
      else if (particle_datasets_[index] == "ParticleDisplacements_y") 
        update_particle_displacements_
          (array_double,nx,ny,nz,
           particle,it,ia,lower_block[1],h4[IY],1);
      else if (particle_datasets_[index] == "ParticleDisplacements_z") 
        update_particle_displacements_
          (array_double,nx,ny,nz,
           particle,it,ia,lower_block[2],h4[IZ],2);
      
    }
  }  
}

//======================================================================

int EnzoInitialHdf5::read_dataset_
(void ** data, int index,
 Block * block,
 std::string file_name,
 std::string axis_map,
 std::string dataset,
 int *mx, int *my, int *mz,
 int *nx, int *ny, int *nz,
 int *gx, int *gy, int *gz,
 int n4[4], double h4[4],
 int *IX, int *IY, int *IZ)
 
{

  // Get the grid size at level_
  Hierarchy * hierarchy = cello::simulation()->hierarchy();
  double lower_domain[3];
  double upper_domain[3];
  hierarchy->lower(lower_domain, lower_domain+1, lower_domain+2);
  hierarchy->upper(upper_domain, upper_domain+1, upper_domain+2);

  double lower_block[3];
  double upper_block[3];
  block->lower(lower_block, lower_block+1, lower_block+2);
  block->upper(upper_block, upper_block+1, upper_block+2);

  Field field = block->data()->field();
  
  field.dimensions (0,mx,my,mz);
  field.size         (nx,ny,nz);
  field.ghost_depth(0,gx,gy,gz);

  // Open the file

  FileHdf5 * file = new FileHdf5 ("./",file_name);
  file->file_open();

  // Read the domain dimensions

  *IX = axis_map.find ("x");
  *IY = axis_map.find ("y");
  *IZ = axis_map.find ("z");

  ASSERT3 ("EnzoInitialHdf5::enforce_block()",
           "bad field coordinates %d %d %d",
           (*IX),(*IY),(*IZ),
           (((*IX)<4)&&((*IY)<4)&&((*IZ)<4)) &&
           (((*IX) != (*IY)) || ((*IY)==-1 && (*IZ) == -1)) &&
           (((*IX) != (*IY) && (*IY) != (*IZ)) || ((*IZ) == -1)));
    
  int m4[4] = {0,0,0,0};
  int type_data = type_unknown;
  file-> data_open (dataset, &type_data,
                    m4,m4+1,m4+2,m4+3);
  // compute cell widths
  for (int i=0; i<4; i++) { h4[i] = 1.0; };
  h4[(*IX)] = (upper_block[0] - lower_block[0]) / (*nx);
  h4[(*IY)] = (upper_block[1] - lower_block[1]) / (*ny);
  h4[(*IZ)] = (upper_block[2] - lower_block[2]) / (*nz);

  // determine offsets
  int o4[4] = {0,0,0,0};
  o4[(*IX)] = (lower_block[0] - lower_domain[0]) / h4[(*IX)];
  o4[(*IY)] = (lower_block[1] - lower_domain[1]) / h4[(*IY)];
  o4[(*IZ)] = (lower_block[2] - lower_domain[2]) / h4[(*IZ)];

  // adjust offsets if domain is larger than file input
  // (e.g. to allow N=1024^3 using 512^3 input files
  // for scaling tests)

  if (o4[(*IX)] >= m4[(*IX)]) o4[(*IX)] = o4[(*IX)] % m4[(*IX)];
  if (o4[(*IY)] >= m4[(*IY)]) o4[(*IY)] = o4[(*IY)] % m4[(*IY)];
  if (o4[(*IZ)] >= m4[(*IZ)]) o4[(*IZ)] = o4[(*IZ)] % m4[(*IZ)];

  for (int i=0; i<4; i++) n4[i] = 1;
  n4[(*IX)] = (upper_block[0] - lower_block[0]) / h4[(*IX)];
  n4[(*IY)] = (upper_block[1] - lower_block[1]) / h4[(*IY)];
  n4[(*IZ)] = (upper_block[2] - lower_block[2]) / h4[(*IZ)];

  // open the dataspace
  file-> data_slice
    (m4[0],m4[1],m4[2],m4[3],
     n4[0],n4[1],n4[2],n4[3],
     o4[0],o4[1],o4[2],o4[3]);

  // create memory space
  // (fields was n4[(*IX)],n4[(*IY)],n4[(*IZ)])
  file->mem_create ((*nx),(*ny),(*nz),(*nx),(*ny),(*nz),0,0,0);
    
  // input domain size

  const int n = (*nx)*(*ny)*(*nz);
  if (type_data == type_single) {
    (*data) = new float [n];
  } else if (type_data == type_double) {
    (*data) = new double [n];
  } else {
    ERROR2 ("EnzoInitialHdf5::enforce_block()",
            "Unsupported data type %d in file %s",
            type_data,file_name.c_str());
  }

  file->data_read ((*data));
  file->data_close();
  file->file_close();
  delete file;

  return type_data;
}

//----------------------------------------------------------------------

template <class T>
void EnzoInitialHdf5::update_particle_displacements_
( T * array,int nx, int ny, int nz,
  Particle particle, int it, int ia, double lower, double h,
  int axis)
{
  if (axis == 0) {
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
        for (int ix=0; ix<nx; ix++) {
          int ip = ix + nx*(iy + ny*iz);
          int ib,io;
          particle.index(ip,&ib,&io);
          array = ( T *) particle.attribute_array(it,ia,ib);
          array[io] += lower + (ix+0.5)*h;
        }
      }
    }
  } else if (axis == 1) {
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
        for (int ix=0; ix<nx; ix++) {
          int ip = ix + nx*(iy + ny*iz);
          int ib,io;
          particle.index(ip,&ib,&io);
          array = ( T *) particle.attribute_array(it,ia,ib);
          array[io] += lower + (iy+0.5)*h;
        }
      }
    }
  } else if (axis == 2) {
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
        for (int ix=0; ix<nx; ix++) {
          int ip = ix + nx*(iy + ny*iz);
          int ib,io;
          particle.index(ip,&ib,&io);
          array = ( T *) particle.attribute_array(it,ia,ib);
          array[io] += lower + (iz+0.5)*h;
        }
      }
    }
  }    
}

//----------------------------------------------------------------------

template <class T>
void EnzoInitialHdf5::copy_field_data_to_array_
(enzo_float * array, T * data,
 int mx,int my,int mz,int nx,int ny,int nz,int gx,int gy,int gz,int n4[4],
 int IX, int IY) const
{
  for (int iz=0; iz<nz; iz++) {
    int jz = iz+gz;
    for (int iy=0; iy<ny; iy++) {
      int jy = iy+gy;
      for (int ix=0; ix<nx; ix++) {
	int jx = ix+gx;
	int i = ix+n4[IX]*(iy+n4[IY]*iz);
	int j = jx+mx*(jy+my*jz);
	array[j] = data[i];
      }
    }
  }
}

//----------------------------------------------------------------------

template <class T, class S>
void EnzoInitialHdf5::copy_particle_data_to_array_
(T * array, S * data,
 Particle particle, int it, int ia, int np)
{
  for (int ip=0; ip<np; ip++) {
    int ib,io;
    particle.index(ip,&ib,&io);
    array = (T*)particle.attribute_array(it,ia,ib);
    array[io] = data[ip];
  }
}
