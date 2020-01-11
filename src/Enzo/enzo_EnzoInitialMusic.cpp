// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitiaHdf5.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-06-23
/// @brief    Read initial conditions from HDF5
///           (multi-scale cosmological initial conditions)
#include "enzo.hpp"
#include <chrono>
#include <thread>

// #define DEBUG_THROTTLE

//----------------------------------------------------------------------

static CmiNodeLock throttle_node_lock;

//----------------------------------------------------------------------
void mutex_init()
{
  throttle_node_lock = CmiCreateLock();
}

//----------------------------------------------------------------------

EnzoInitialMusic::EnzoInitialMusic
(int cycle,
 double time,
 const EnzoConfig * enzo_config,
 int level) throw()
  : Initial(cycle,time),
    level_(level),
    field_files_     (enzo_config->initial_music_field_files),
    field_datasets_  (enzo_config->initial_music_field_datasets),
    field_coords_    (enzo_config->initial_music_field_coords),
    field_names_     (enzo_config->initial_music_field_names),
    particle_files_     (enzo_config->initial_music_particle_files),
    particle_datasets_  (enzo_config->initial_music_particle_datasets),
    particle_coords_    (enzo_config->initial_music_particle_coords),
    particle_types_     (enzo_config->initial_music_particle_types),
    particle_attributes_(enzo_config->initial_music_particle_attributes),
    throttle_internode_ (enzo_config->initial_music_throttle_internode),
    throttle_intranode_ (enzo_config->initial_music_throttle_intranode),
    throttle_group_size_    (enzo_config->initial_music_throttle_group_size),
    throttle_seconds_stagger_ (enzo_config->initial_music_throttle_seconds_stagger),
    throttle_seconds_delay_ (enzo_config->initial_music_throttle_seconds_delay)
{
#ifndef CONFIG_SMP_MODE
  ASSERT ("EnzoInitialMusic::EnzoInitialMusic",
          "CONFIG_SMP_MODE must be enabled for intranode throttling",
          (! throttle_intranode_));
#endif          
}

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

  p | throttle_intranode_;
  p | throttle_internode_;
  p | throttle_group_size_;
  p | throttle_seconds_stagger_;
  p | throttle_seconds_delay_;
}

//----------------------------------------------------------------------

void EnzoInitialMusic::enforce_block
( Block * block, const Hierarchy * hierarchy ) throw()
{

  if (block->level() != level_) return;

  // Optionally pause before reading if throttling enabled.  For
  // reducing filesystem contention on large runs
  throttle_stagger_();
 
  // Get the grid size at level_
  double lower_domain[3];
  double upper_domain[3];
  hierarchy->lower(lower_domain, lower_domain+1, lower_domain+2);
  hierarchy->upper(upper_domain, upper_domain+1, upper_domain+2);

  double lower_block[3];
  double upper_block[3];
  block->lower(lower_block, lower_block+1, lower_block+2);
  block->upper(upper_block, upper_block+1, upper_block+2);

  Field field = block->data()->field();

  for (size_t index=0; index<field_files_.size(); index++) {

    std::string file_name = field_files_[index];

    // Block size

    int mx,my,mz;
    int nx,ny,nz;
    int gx,gy,gz;
 
    field.dimensions (0,&mx,&my,&mz);
    field.size         (&nx,&ny,&nz);
    field.ghost_depth(0,&gx,&gy,&gz);

    // Open the field file

    FileHdf5 file("./",file_name);

    if (throttle_intranode_) {
      CmiLock(throttle_node_lock);
    }

#ifdef DEBUG_THROTTLE
    CkPrintf ("%d %g DEBUG_THROTTLE opening %s\n",
              CkMyPe(),cello::simulation()->timer(),file_name.c_str());
    fflush(stdout);
#endif
      
    file.file_open();

    // Read the domain dimensions

    const int IX = field_coords_[index].find ("x");
    const int IY = field_coords_[index].find ("y");
    const int IZ = field_coords_[index].find ("z");

    ASSERT3 ("EnzoInitialMusic::enforce_block()",
	    "bad field coordinates %d %d %d",
	     IX,IY,IZ,
	     ((IX<4)&&(IY<4)&&(IZ<4)) &&
	     ((IX != IY) || (IY==-1 && IZ == -1)) &&
	     ((IX != IY && IY != IZ) || (IZ == -1)));
    
    int m4[4] = {0};
    int type_data = type_unknown;
    file. data_open (field_datasets_[index], &type_data,
		     m4,m4+1,m4+2,m4+3);
    // compute cell widths
    double h4[4] = {1};
    h4[IX] = (upper_block[0] - lower_block[0]) / nx;
    h4[IY] = (upper_block[1] - lower_block[1]) / ny;
    h4[IZ] = (upper_block[2] - lower_block[2]) / nz;

    // determine offsets
    int o4[4] = {0};
    o4[IX] = (lower_block[0] - lower_domain[0]) / h4[IX];
    o4[IY] = (lower_block[1] - lower_domain[1]) / h4[IY];
    o4[IZ] = (lower_block[2] - lower_domain[2]) / h4[IZ];

    // adjust offsets if domain is larger than file input
    // (e.g. to allow N=1024^3 using 512^3 input files
    // for scaling tests)

    if (o4[IX] >= m4[IX]) o4[IX] = o4[IX] % m4[IX];
    if (o4[IY] >= m4[IY]) o4[IY] = o4[IY] % m4[IY];
    if (o4[IZ] >= m4[IZ]) o4[IZ] = o4[IZ] % m4[IZ];

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
    union {
      void   * data;
      float  * data_float;
      double * data_double;
    };
    
    if (type_data == type_single) {
      data_float = new float [nx*ny*nz];
    } else if (type_data == type_double) {
      data_double = new double [nx*ny*nz];
    } else {
      ERROR3 ("EnzoInitialMusic::enforce_block()",
	      "Unsupported data type %d in file %s field %s",
	      type_data,file_name.c_str(),field_datasets_[index].c_str());
    }

    file.data_read (data);

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
    
    file.data_close();
    file.file_close();
#ifdef DEBUG_THROTTLE
      CkPrintf ("%d %g DEBUG_THROTTLE closed %s\n",
                CkMyPe(),cello::simulation()->timer(),file_name.c_str());
      fflush(stdout);
#endif    
      throttle_delay_();
      if (throttle_intranode_) {
	CmiUnlock(throttle_node_lock);
      }


  }

  for (size_t index=0; index<particle_files_.size(); index++) {

    std::string file_name = particle_files_[index];

    // Open the particle file

    if (throttle_intranode_) {
      CmiLock(throttle_node_lock);
    }
#ifdef DEBUG_THROTTLE
    CkPrintf ("%d %g DEBUG_THROTTLE opening %s\n",
              CkMyPe(),cello::simulation()->timer(),file_name.c_str());
    fflush(stdout);
#endif


    FileHdf5 file("./",file_name);

    file.file_open();

    // Open the dataset
    int m4[4] = {0};
    int type_data = type_unknown;
    file. data_open (particle_datasets_[index], &type_data,
		     m4,m4+1,m4+2,m4+3);

    // Block size

    int mx,my,mz;
    int nx,ny,nz;
    int gx,gy,gz;
 
    field.dimensions (0,&mx,&my,&mz);
    field.size         (&nx,&ny,&nz);
    field.ghost_depth(0,&gx,&gy,&gz);

    // coordinate mapping
    const int IX = particle_coords_[index].find ("x");
    const int IY = particle_coords_[index].find ("y");
    const int IZ = particle_coords_[index].find ("z");

    // compute cell widths
    double h4[4] = {1};
    h4[IX] = (upper_block[0] - lower_block[0]) / nx;
    h4[IY] = (upper_block[1] - lower_block[1]) / ny;
    h4[IZ] = (upper_block[2] - lower_block[2]) / nz;

    // determine offsets
    int o4[4] = {0};
    o4[IX] = (lower_block[0] - lower_domain[0]) / h4[IX];
    o4[IY] = (lower_block[1] - lower_domain[1]) / h4[IY];
    o4[IZ] = (lower_block[2] - lower_domain[2]) / h4[IZ];

    if (o4[IX] >= m4[IX]) o4[IX] = o4[IX] % m4[IX];
    if (o4[IY] >= m4[IY]) o4[IY] = o4[IY] % m4[IY];
    if (o4[IZ] >= m4[IZ]) o4[IZ] = o4[IZ] % m4[IZ];

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

    file.mem_create (nx,ny,nz,nx,ny,nz,0,0,0);

    // input domain size
    union {
      void   * data;
      float  * data_float;
      double * data_double;
    };

    if (type_data == type_single) {
      data_float = new float [nx*ny*nz];
    } else if (type_data == type_double) {
      data_double = new double [nx*ny*nz];
    } else {
      ERROR3 ("EnzoInitialMusic::enforce_block()",
	      "Unsupported data type %d in file %s particle dataset %s",
	      type_data,file_name.c_str(),particle_datasets_[index].c_str());
    }
    
    // read data and close file
    file.data_read (data);
    file.data_close();
    file.file_close();

#ifdef DEBUG_THROTTLE
    CkPrintf ("%d %g DEBUG_THROTTLE closed %s\n",
              CkMyPe(),cello::simulation()->timer(),file_name.c_str());
    fflush(stdout);
#endif    

    throttle_delay_();
    if (throttle_intranode_) {
      CmiUnlock(throttle_node_lock);
    }
    
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
      ERROR3 ("EnzoInitialMusic::enforce_block()",
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
      
      if (particle_datasets_[index] == "ParticleDisplacements_x") {
	for (int iz=0; iz<nz; iz++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int ix=0; ix<nx; ix++) {
	      int ip = ix + nx*(iy + ny*iz);
	      int ib,io;
	      particle.index(ip,&ib,&io);
	      array = particle.attribute_array(it,ia,ib);
	      array_float[io] += lower_block[0] + (ix+0.5)*h4[IX];
	    }
	  }
	}
      } else if (particle_datasets_[index] == "ParticleDisplacements_y") {
	for (int iz=0; iz<nz; iz++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int ix=0; ix<nx; ix++) {
	      int ip = ix + nx*(iy + ny*iz);
	      int ib,io;
	      particle.index(ip,&ib,&io);
	      array = particle.attribute_array(it,ia,ib);
	      array_float[io] += lower_block[1] + (iy+0.5)*h4[IY];
	    }
	  }
	}
      } else if (particle_datasets_[index] == "ParticleDisplacements_z") {
	for (int iz=0; iz<nz; iz++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int ix=0; ix<nx; ix++) {
	      int ip = ix + nx*(iy + ny*iz);
	      int ib,io;
	      particle.index(ip,&ib,&io);
	      array = particle.attribute_array(it,ia,ib);
	      array_float[io] += lower_block[2] + (iz+0.5)*h4[IZ];
	    }
	  }
	}
      }

    } else { // (type_array != type_single) {
      
      if (particle_datasets_[index] == "ParticleDisplacements_x") {
	for (int iz=0; iz<nz; iz++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int ix=0; ix<nx; ix++) {
	      int ip = ix + nx*(iy + ny*iz);
	      int ib,io;
	      particle.index(ip,&ib,&io);
	      array = particle.attribute_array(it,ia,ib);
	      array_double[io] += lower_block[0] + (ix+0.5)*h4[IX];
	    }
	  }
	}
      } else if (particle_datasets_[index] == "ParticleDisplacements_y") {
	for (int iz=0; iz<nz; iz++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int ix=0; ix<nx; ix++) {
	      int ip = ix + nx*(iy + ny*iz);
	      int ib,io;
	      particle.index(ip,&ib,&io);
	      array = particle.attribute_array(it,ia,ib);
	      array_double[io] += lower_block[1] + (iy+0.5)*h4[IY];
	    }
	  }
	}
      } else if (particle_datasets_[index] == "ParticleDisplacements_z") {
	for (int iz=0; iz<nz; iz++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int ix=0; ix<nx; ix++) {
	      int ip = ix + nx*(iy + ny*iz);
	      int ib,io;
	      particle.index(ip,&ib,&io);
	      array = particle.attribute_array(it,ia,ib);
	      array_double[io] += lower_block[2] + (iz+0.5)*h4[IZ];
	    }
	  }
	}
      }
    }
  }  
}

//----------------------------------------------------------------------

void EnzoInitialMusic::throttle_stagger_()
{
  if (throttle_internode_) {
    //--------------------------------------------------
    const int in= CkMyPe()/CONFIG_NODE_SIZE;
    static bool not_first[CONFIG_NODE_SIZE] = {false};
    if (throttle_group_size_ != 0 && ! not_first[in]) {
      not_first[in] = true;
      int ms = 1000*((CkMyPe() % throttle_group_size_) * throttle_seconds_stagger_);
#ifdef DEBUG_THROTTLE  
      CkPrintf ("DEBUG_THROTTLE %d %g %d ms stagger\n",CkMyPe(),cello::simulation()->timer(),ms);
      fflush(stdout);
#endif      
      std::this_thread::sleep_for(std::chrono::milliseconds(ms));
    }
  }
}

//----------------------------------------------------------------------

void EnzoInitialMusic::throttle_delay_()
{
  if (throttle_internode_) {
    int ms = 1000*throttle_seconds_delay_;
#ifdef DEBUG_THROTTLE  
    CkPrintf ("DEBUG_THROTTLE %d %g %d ms delay\n",
	      CkMyPe(),cello::simulation()->timer(),ms);
    fflush(stdout);
#endif      
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
  }
}

//======================================================================

template <class T>
void EnzoInitialMusic::copy_field_data_to_array_
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
void EnzoInitialMusic::copy_particle_data_to_array_
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
