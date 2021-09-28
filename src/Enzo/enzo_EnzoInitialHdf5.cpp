/// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitiaHdf5.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-05-13
/// @brief    Read initial conditions from HDF5
#include "enzo.hpp"
#include <chrono>
#include <thread>

#define CHECK_COSMO_PARAMS true

//----------------------------------------------------------------------

EnzoInitialHdf5::EnzoInitialHdf5
(int cycle,
 double time,
 int max_level,
 std::string                 format,
 const int                   blocking[3],
 std::vector < std::string > field_files,
 std::vector < std::string > field_datasets,
 std::vector < std::string > field_coords,
 std::vector < std::string > field_names,
 std::vector < std::string > particle_files,
 std::vector < std::string > particle_datasets,
 std::vector < std::string > particle_coords,
 std::vector < std::string > particle_types,
 std::vector < std::string > particle_attributes) throw()
   : Initial (cycle,time),
     max_level_(max_level),
     format_ (format),
     field_files_ (field_files),
     field_datasets_ (field_datasets),
     field_coords_ (field_coords),
     field_names_ (field_names),
     particle_files_ (particle_files),
     particle_datasets_ (particle_datasets),
     particle_coords_ (particle_coords),
     particle_types_ (particle_types),
     particle_attributes_ (particle_attributes),
     l_particle_displacements_(false),
     particle_position_names_()
{
  for (int i=0; i<3; i++) blocking_[i]=blocking[i];

  if (format == "music") {
    l_particle_displacements_ = true;
    particle_position_names_[0] = "ParticleDisplacements_x";
    particle_position_names_[1] = "ParticleDisplacements_y";
    particle_position_names_[2] = "ParticleDisplacements_z";
    for (int i=0; i<field_coords.size(); i++)  field_coords_[i] = "tzyx";
    for (int i=0; i<particle_coords.size(); i++) particle_coords_[i] = "tzyx";

  } else {
    ERROR1 ("EnzoInitialHdf5::EnzoInitialHdf5()",
            "Unsupported format '%s'",
            format.c_str());
  }
  i_sync_msg_ = cello::scalar_descr_sync()->new_value("initial_hdf5:msg");
}

//----------------------------------------------------------------------

void EnzoInitialHdf5::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change

  p | max_level_;

  p | field_files_;
  p | field_datasets_;
  p | field_coords_;
  p | field_names_;

  p | particle_files_;
  p | particle_datasets_;
  p | particle_coords_;
  p | particle_types_;
  p | particle_attributes_;
  PUParray (p,particle_position_names_,3);
  p | l_particle_displacements_;

}

//----------------------------------------------------------------------

void EnzoInitialHdf5::enforce_block
( Block * block, const Hierarchy * hierarchy_unused ) throw()
{
  if (! (0 <= block->level() && block->level() <= max_level_) ) {
    // if level not in range, then return and call initial_done()
    block->initial_done();
    return;
  } else if (! is_reader_(block->index())) {
    // else if not reader, will be expected to receive data from a reader, so
    // return and sit back and wait, but /don't/ call initial_done() until
    // your reader says you're ready
    Sync * sync_msg = psync_msg_(block);
    return;
  }

  // Assert: to reach this point, block must be a reading block

  int array_lower[3],array_upper[3];
  root_block_range_(block->index(),array_lower,array_upper);

  Field field = block->data()->field();

  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;
  int n4[4],IX,IY,IZ;
  double h4[4];

  field.dimensions (0,&mx,&my,&mz);
  field.size         (&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);

  double lower_block[3],upper_block[3];
  block->lower(lower_block,lower_block+1,lower_block+2);
  block->upper(upper_block,upper_block+1,upper_block+2);
  static std::map<std::string,int> close_count;

  // Maintain running count of messages sent
  int count_messages = 0;

  // Read in Field files
  for (size_t index=0; index<field_files_.size(); index++) {

    // Open the file

    FileHdf5 * file = new FileHdf5 ("./",field_files_[index]);
    file->file_open();

    // Double-check cosmology parameters if CHECK_COSMO_PARAMS is true
    if (CHECK_COSMO_PARAMS) {
      check_cosmology_(file);
    }

    // Open the dataset

    int m4[4] = {0,0,0,0};
    int type_data = type_unknown;
    file-> data_open (field_datasets_[index], &type_data,
                      m4,m4+1,m4+2,m4+3);

    ASSERT1("EnzoInitialHdf5::enforce_block()",
           "Unsupported type_data %d",
            type_data,
            ( (type_data == type_single) ||
              (type_data == type_double) ) );

    // Count number of messages to send per block
    ++count_messages;

    // Loop over root-level blocks in range of this reader
    for (int ax=array_lower[0]; ax<array_upper[0]; ax++) {
      for (int ay=array_lower[1]; ay<array_upper[1]; ay++) {
        for (int az=array_lower[2]; az<array_upper[2]; az++) {
          int block_index[3] = {ax,ay,az};
          Index index_block(ax,ay,az);

          char * data;
          read_dataset_
            (file, &data,index_block, type_data,
             lower_block,upper_block,block_index,
             field_coords_[index],
             nx,ny,nz,m4,n4,h4,&IX,&IY,&IZ);

          if (index_block == block->index() ) {

            // local block: copy directly to field
            copy_dataset_to_field_
              (block, field_names_[index],type_data,
               data,mx,my,mz,nx,ny,nz,gx,gy,gz,n4,IX,IY);

          } else {

            // remote block: pack message and send
            MsgInitial * msg_initial = new MsgInitial;
            msg_initial->set_dataset (n4,h4,nx,ny,nz,IX,IY,IZ);
            msg_initial->set_field_data
              (field_names_[index],data,nx*ny*nz,type_data);
            enzo::block_array()[index_block].p_initial_hdf5_recv(msg_initial);
          }

          delete_array_ (&data,type_data);

        }
      }
    }
    file->data_close();
    file->file_close();
    delete file;
  }

  // Read in particle files

  for (size_t index=0; index<particle_files_.size(); index++) {

    FileHdf5 * file = new FileHdf5 ("./",particle_files_[index]);
    file->file_open();

    // Double-check cosmology parameters if CHECK_COSMO_PARAMS is true
    if (CHECK_COSMO_PARAMS) {
      check_cosmology_(file);
    }

    // Open the dataset

    int m4[4] = {0,0,0,0};
    int type_data = type_unknown;
    file-> data_open (particle_datasets_[index], &type_data,
                      m4,m4+1,m4+2,m4+3);

    ASSERT1("EnzoInitialHdf5::enforce_block()",
           "Unsupported type_data %d",
            type_data,
            ( (type_data == type_single) ||
              (type_data == type_double) ) );

    // Count number of messages to send per block
    ++count_messages;

    // Loop over root-level blocks in range of this reader
    for (int ax=array_lower[0]; ax<array_upper[0]; ax++) {
      for (int ay=array_lower[1]; ay<array_upper[1]; ay++) {
        for (int az=array_lower[2]; az<array_upper[2]; az++) {

          int block_index[3] = {ax,ay,az};
          Index index_block(ax,ay,az);

          char * data;
          read_dataset_
            (file, &data,index_block, type_data,
             lower_block,upper_block,block_index,
             particle_coords_[index],
             nx,ny,nz,m4,n4,h4,&IX,&IY,&IZ);

          if (index_block == block->index() ) {

            // local block: copy directly to particle
            copy_dataset_to_particle_
              (block,
               particle_types_[index],
               particle_attributes_[index],
               type_data,
               data,
               nx,ny,nz,
               h4,IX,IY,IZ);

          } else {

            // remote block: pack message and send
            MsgInitial * msg_initial = new MsgInitial;
            msg_initial->set_dataset (n4,h4,nx,ny,nz,IX,IY,IZ);
            msg_initial->set_particle_data
              (particle_types_[index],
               particle_attributes_[index],
               data,nx*ny*nz,type_data);

            enzo::block_array()[index_block].p_initial_hdf5_recv(msg_initial);

          }
          delete_array_ (&data,type_data);
        }
      }
    }
    file->data_close();
    file->file_close();
    delete file;
  }

  for (int ax=array_lower[0]; ax<array_upper[0]; ax++) {
    for (int ay=array_lower[1]; ay<array_upper[1]; ay++) {
      for (int az=array_lower[2]; az<array_upper[2]; az++) {
        Index index_block(ax,ay,az);
        if (index_block != block->index() ) {
          MsgInitial * msg_initial = new MsgInitial;
          msg_initial->set_count(count_messages + 1);
          // send empty message with count of number of messages sent
          // (including this one)
          enzo::block_array()[index_block].p_initial_hdf5_recv(msg_initial);
        }
      }
    }
  }
  block->initial_done();
}

//----------------------------------------------------------------------

void EnzoBlock::p_initial_hdf5_recv(MsgInitial * msg_initial)
{
  EnzoInitialHdf5 * initial = static_cast<EnzoInitialHdf5*> (this->initial());
  initial->recv_data(this,msg_initial);
}

//----------------------------------------------------------------------

void EnzoInitialHdf5::recv_data (Block * block, MsgInitial * msg_initial)
{
  // Exit when count reached (set_stop() may be called at any time)
  Sync * sync_msg = psync_msg_(block);
  int count = msg_initial->count();
  if ( count > 0) {
    sync_msg->set_stop(count);
  }

  // Copy data from message to block data
  if (msg_initial->data_type() == "field") {
    // extract parameters from MsgInitial
    Field field = block->data()->field();
    int n4[4];
    double h4[4];
    int nx,ny,nz;
    int IX,IY,IZ;
    std::string field_name;
    char * data;
    int data_precision;
    msg_initial->get_dataset (n4,h4,&nx,&ny,&nz,&IX,&IY,&IZ);
    msg_initial->get_field_data
              (&field_name,&data,&data_precision);
    int index_field = field.field_id(field_name);
    int mx,my,mz;
    int gx,gy,gz;
    field.dimensions (index_field,&mx,&my,&mz);
    field.ghost_depth(0,&gx,&gy,&gz);

    // Copy data to field
    copy_dataset_to_field_
      (block, field_name,data_precision,
       data,mx,my,mz,nx,ny,nz,gx,gy,gz,n4,IX,IY);

  } else if (msg_initial->data_type() == "particle") {

    // extract parameters from MsgInitial
    Particle particle = block->data()->particle();
    int n4[4];
    double h4[4];
    int nx,ny,nz;
    int IX,IY,IZ;
    msg_initial->get_dataset (n4,h4,&nx,&ny,&nz,&IX,&IY,&IZ);

    std::string particle_type, particle_attribute;
    char * data;
    int data_size, data_precision;
    msg_initial->get_particle_data
              (&particle_type,
               &particle_attribute,
               &data,&data_size,&data_precision);

    // Copy data to particle
    copy_dataset_to_particle_
      (block,
       particle_type,
       particle_attribute,
       data_precision,
       (char *)data,
       nx,ny,nz,
       h4,IX,IY,IZ);

  }

  if (sync_msg->next()) {
    // reset for next call (note not resetting at start since may get
    // called after messages received)
    sync_msg->reset();
    block->initial_done();
  }
  delete msg_initial;
}

//======================================================================

void EnzoInitialHdf5::copy_dataset_to_field_
(Block * block,
 std::string field_name, int type_data,
 char * data,
 int mx, int my, int mz,
 int nx, int ny, int nz,
 int gx, int gy, int gz,
 int n4[4], int IX, int IY)
{
  Field field = block->data()->field();

  // Destination is this block--copy directly
  enzo_float * array = (enzo_float *) field.values(field_name);

  if (type_data == type_single) {
    copy_field_data_to_array_
      (array,(float *)data,mx,my,mz,nx,ny,nz,gx,gy,gz,n4,IX,IY);
  } else if (type_data == type_double) {
    copy_field_data_to_array_
      (array,(double *)data,mx,my,mz,nx,ny,nz,gx,gy,gz,n4,IX,IY);
  }
}

//----------------------------------------------------------------------

void EnzoInitialHdf5::copy_dataset_to_particle_
(Block * block,
 std::string particle_type,
 std::string particle_attribute,
 int type_data,
 char * data,
 int nx, int ny, int nz,
 double h4[4], int IX, int IY, int IZ)
{

  // Create particles and initialize them
  Particle particle = block->data()->particle();

  const int it = particle.type_index(particle_type);
  const int ia = particle.attribute_index(it,particle_attribute);

  const int np = nx*ny*nz;

  // insert particles if they don't exist yet
  // (assumes dataset covers entire block)
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
        (array_float,(float *)data,particle,it,ia,np);
    } else if (type_data == type_double) {
      copy_particle_data_to_array_
        (array_float,(double *)data,particle,it,ia,np);
    }
  } else if (type_array == type_double) {
    if (type_data == type_single) {
      copy_particle_data_to_array_
        (array_double,(float *)data,particle,it,ia,np);
    } else if (type_data == type_double) {
      copy_particle_data_to_array_
        (array_double,(double *)data,particle,it,ia,np);
    }
  } else {
    ERROR3 ("EnzoInitialHdf5::enforce_block()",
            "Unsupported particle precision %s for "
            "particle type %s attribute %s",
            cello::precision_name[type_array],
            particle.type_name(it).c_str(),
            particle.attribute_name(it,ia).c_str());
  }

  // update positions with displacements if needed

  double lower_block[3];
  block->lower(lower_block,lower_block+1,lower_block+2);

  if (l_particle_displacements_) {
    for (int axis=0; axis<3; axis++) {
      // if attribute corresponds to position along the given axis,
      // add to the attribute values the offsets of the cell centers
      // along that axis
      if (ia == particle.attribute_position(it,axis)) {
        if (type_array == type_single) {
          update_particle_displacements_
            (array_float,nx,ny,nz,
             particle,it,ia,lower_block[axis],h4[IX],axis);
        } else {
          update_particle_displacements_
            (array_double, nx,ny,nz,
             particle,it,ia,lower_block[axis],h4[IX],axis);
        }
      }
    }
  }
}

//----------------------------------------------------------------------

int EnzoInitialHdf5::is_reader_ (Index index)
{
  int a3[3];
  index.array(a3,a3+1,a3+2);
  const int level = index.level();
  return ( (level == 0) &&
           ( a3[0] % blocking_[0] == 0) &&
           ( a3[1] % blocking_[1] == 0) &&
           ( a3[2] % blocking_[2] == 0));
}

//----------------------------------------------------------------------

void EnzoInitialHdf5::root_block_range_(Index index, int array_lower[3], int array_upper[3])
{
  // Get array-of-octrees blocking
  int root_blocks[3];
  cello::hierarchy()->root_blocks
    (root_blocks,  root_blocks+1, root_blocks+2);
  // Get this (reader) block index in root blocking
  index.array(array_lower, array_lower+1,array_lower+2);
  array_upper[0] = std::min(array_lower[0] + blocking_[0],root_blocks[0]);
  array_upper[1] = std::min(array_lower[1] + blocking_[1],root_blocks[1]);
  array_upper[2] = std::min(array_lower[2] + blocking_[2],root_blocks[2]);
}

//----------------------------------------------------------------------

void EnzoInitialHdf5::read_dataset_
(File * file, char ** data, Index index, int type_data,
 double lower_block[3], double upper_block[3],
 int block_index[3],
 std::string axis_map,
 int nx, int ny, int nz,
 int m4[4], int n4[4],double h4[4],
 int *IX, int *IY, int *IZ)
{

  // Get the grid size at level_
  Hierarchy * hierarchy = cello::simulation()->hierarchy();
  double lower_domain[3];
  double upper_domain[3];
  hierarchy->lower(lower_domain, lower_domain+1, lower_domain+2);
  hierarchy->upper(upper_domain, upper_domain+1, upper_domain+2);

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

  // field size
  n4[0] = n4[1] = n4[2] = n4[3] = 1;
  n4[(*IX)] = nx;
  n4[(*IY)] = ny;
  n4[(*IZ)] = nz;

  // compute cell widths
  h4[0] = h4[1] = h4[2] = h4[3] = 1.0;
  h4[(*IX)] = (upper_block[0] - lower_block[0]) / nx;
  h4[(*IY)] = (upper_block[1] - lower_block[1]) / ny;
  h4[(*IZ)] = (upper_block[2] - lower_block[2]) / nz;

  // determine offsets
  int o4[4] = {0,0,0,0};
  o4[(*IX)] = block_index[0]*nx;
  o4[(*IY)] = block_index[1]*ny;
  o4[(*IZ)] = block_index[2]*nz;

  // open the dataspace
  file-> data_slice
    (m4[0],m4[1],m4[2],m4[3],
     n4[0],n4[1],n4[2],n4[3],
     o4[0],o4[1],o4[2],o4[3]);

  // create memory space
  // (fields was n4[(*IX)],n4[(*IY)],n4[(*IZ)])
  file->mem_create (nx,ny,nz,nx,ny,nz,0,0,0);

  // input domain size

  const int n = nx*ny*nz;
  (*data) = allocate_array_ (n,type_data);

  file->data_read ((*data));
}

//----------------------------------------------------------------------

template <class T>
void EnzoInitialHdf5::update_particle_displacements_
( T * array,int nx, int ny, int nz,
  Particle particle, int it, int ia, double lower, double h,
  int axis)
{
  // set (bx,by,bz) = e[axis] to avoid multiple loops and conditionals
  // inside loops
  const int bx = (axis == 0) ? 1 : 0;
  const int by = (axis == 1) ? 1 : 0;
  const int bz = (axis == 2) ? 1 : 0;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
        int ip = ix + nx*(iy + ny*iz);
        int ib,io;
        // get batch and offset into batch for the ip'th particle
        particle.index(ip,&ib,&io);
        array = ( T *) particle.attribute_array(it,ia,ib);
        array[io] += lower + h*(bx*ix+by*iy+bz*iz + 0.5);
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

//----------------------------------------------------------------------

void EnzoInitialHdf5::check_cosmology_(File * file) const
{
  EnzoPhysicsCosmology * cosmo = enzo::cosmology();

  if (cosmo != nullptr) {
    if (format_ == "music") {

      // Verify cosmology parameters in HDF5 file match those in Enzo-E
      // parameter file to machine precision assumes "MUSIC" format with
      // single-precision attributes To disable, set CHECK_COSMO_PARAMS
      // define to "false" above

      int type = type_unknown;
      float dx,h0,omega_b,omega_m,omega_v,vfact;
      file->file_read_scalar(&dx, "dx", &type);
      file->file_read_scalar(&h0, "h0", &type);
      file->file_read_scalar(&omega_b, "omega_b", &type);
      file->file_read_scalar(&omega_m, "omega_m", &type);
      file->file_read_scalar(&omega_v, "omega_v", &type);
      file->file_read_scalar(&vfact, "vfact", &type);
      double roundoff = cello::machine_epsilon(precision_single);
      ASSERT2("EnzoInitialHdf5",
              "Mismatch in cosmology H0: parameter file %g I.C.'s %g",
              cosmo->hubble_constant_now(),h0,
              cello::err_rel
              (cosmo->hubble_constant_now(),enzo_float(h0)) < roundoff);
      ASSERT2("EnzoInitialHdf5",
              "Mismatch in cosmology omega_b: parameter file %g I.C.'s %g",
              cosmo->omega_baryon_now(),omega_b,
              cello::err_rel
              (cosmo->omega_baryon_now(),enzo_float(omega_b)) < roundoff);
      ASSERT2("EnzoInitialHdf5",
              "Mismatch in cosmology omega_m: parameter file %g I.C.'s %g",
              cosmo->omega_matter_now(),omega_m,
              cello::err_rel
              (cosmo->omega_matter_now(),enzo_float(omega_m)) < roundoff);
      ASSERT2("EnzoInitialHdf5",
              "Mismatch in cosmology omega_v: parameter file %g I.C.'s %g",
              cosmo->omega_lambda_now(),omega_v,
              cello::err_rel
              (cosmo->omega_lambda_now(),enzo_float(omega_v)) < roundoff);
    }
  }
}
