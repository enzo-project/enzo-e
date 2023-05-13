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
 int                         monitor_iter_,
 std::vector < std::string > field_files,
 std::vector < std::string > field_datasets,
 std::vector < std::string > field_coords,
 std::vector < std::string > field_names,
 std::vector < int >         field_levels,
 std::vector < std::string > particle_files,
 std::vector < std::string > particle_datasets,
 std::vector < std::string > particle_coords,
 std::vector < std::string > particle_types,
 std::vector < std::string > particle_attributes,
 std::vector < int >         particle_levels) throw()
   : Initial (cycle,time),
     max_level_(max_level),
     format_ (format),
     monitor_iter_(monitor_iter_),
     field_files_ (field_files),
     field_datasets_ (field_datasets),
     field_coords_ (field_coords),
     field_names_ (field_names),
     field_levels_(field_levels),
     particle_files_ (particle_files),
     particle_datasets_ (particle_datasets),
     particle_coords_ (particle_coords),
     particle_types_ (particle_types),
     particle_attributes_ (particle_attributes),
     particle_levels_(particle_levels),
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
  p | format_;
  PUParray (p,blocking_,3);
  p | monitor_iter_;

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

  FieldLoader field_loader(block, format_);
  ParticleLoader particle_loader(block, l_particle_displacements_, format_);
  int min_level = cello::hierarchy()->min_level();

  // Maintain running count of messages sent
  int count_messages[max_level_ + 1] = {};

  // Read in Field files
  for (size_t index=0; index<field_files_.size(); index++) {
    int level = field_levels_[index];
    if (level > max_level_) {continue;}
    field_loader.open_file(field_files_[index],
                           field_names_[index],
                           field_coords_[index],
                           field_datasets_[index]);
    load_data(count_messages[level],
              block,
              level,
              min_level,
              field_loader);
  }

  // Read in particle files
  for (size_t index=0; index<particle_files_.size(); index++) {
    int level = particle_levels_[index];
    if (level > max_level_) {continue;}
    particle_loader.open_file(particle_files_[index],
                              particle_types_[index],
                              particle_attributes_[index],
                              particle_coords_[index],
                              particle_datasets_[index]);
    load_data(count_messages[level],
              block,
              level,
              min_level,
              particle_loader);
  }

  // Update all blocks in range of this reader with the number of messages sent to them.
  int lower[3], upper[3];
  for (int level = 0; level <= max_level_; level++){
    get_reader_range(block->index(), lower, upper, level);
    for (int ax = lower[0]; ax < upper[0]; ax++) {
      for (int ay = lower[1]; ay < upper[1]; ay++) {
        for (int az = lower[2]; az < upper[2]; az++) {

          Index index_block = block->index_from_global(ax, ay, az, level, min_level);
          if (index_block != block->index()) {
            MsgInitial * msg_initial = new MsgInitial;
            msg_initial->set_count(count_messages[level] + 1);
            enzo::block_array()[index_block].p_initial_hdf5_recv(msg_initial);
          }
        }
      }
    }
  }

  initialize_particle_mass(block);
  block->initial_done();
}

void EnzoInitialHdf5::load_data(int & count_messages,
                                Block * block,
                                int level, int min_level,
                                DataLoader & loader) {
  // Count number of messages to send per block
  ++count_messages;

  int lower[3], upper[3];
  get_reader_range(block->index(), lower, upper, level);
  int region_lower[3];
  cello::hierarchy()->refined_region_lower(region_lower, level-1);

  // Loop over blocks in range of this reader at the given level and use
  // the provided DataLoader to load the associated data onto each of them.
  for (int ax = lower[0]; ax < upper[0]; ax++) {
    for (int ay = lower[1]; ay < upper[1]; ay++) {
      for (int az = lower[2]; az < upper[2]; az++) {
        int block_index[3] = {ax, ay, az};
        for (int i = 0; i < 3; i++) block_index[i] -= (region_lower[i] << 1);
        Index index_block = block->index_from_global(ax, ay, az, level, min_level);
        loader.load(block_index, index_block);
      }
    }
  }
  loader.close_file();
}

void EnzoInitialHdf5::get_reader_range(Index reader_index, int lower[3], int upper[3], int level) throw() {
  // Get the lower and upper values defining the region of parent blocks at level-1
  // which refine to create the blocks on the given level.
  int region_lower[3], region_upper[3];
  cello::hierarchy()->refined_region_lower(region_lower, level-1);
  cello::hierarchy()->refined_region_upper(region_upper, level-1);

  // Get the root level range of the reader block whose index is `reader_index`.
  int root_range_lower[3], root_range_upper[3];
  root_block_range_(reader_index, root_range_lower, root_range_upper);

  // Compute the range of the specfied reader block at the given level.
  for (int i=0; i < 3; i++)
    lower[i] = std::max(root_range_lower[i] << level, region_lower[i] << 1);
  for (int i=0; i < 3; i++)
    upper[i] = std::min(root_range_upper[i] << level, region_upper[i] << 1);
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

  /// Monitor input progress if monitor_iter_ != 0
  static int count_monitor = 0;
  static int count_monitor_out = 0;
  const int blocking = (blocking_[0]*blocking_[1]*blocking_[2]-1);
  if (monitor_iter_ &&
      (msg_initial->data_type()!="field" &&
       msg_initial->data_type()!="particle") &&
      ((count_monitor == 0 || count_monitor == count-1) ||
       ((count_monitor % (monitor_iter_*blocking)) == 0))) {
    cello::monitor()->print("Initial", "hdf5 %d / %d",
                            count_monitor_out,blocking);
    count_monitor_out++;
  }
  count_monitor++;
  
  // Copy data from message to block data
  if (msg_initial->data_type() == "field") {
    FieldLoader field_loader(block, format_);
    char * data;
    field_loader.read_msg(msg_initial, &data);
    field_loader.copy_data_local(data);

  } else if (msg_initial->data_type() == "particle") {
    ParticleLoader particle_loader(block, l_particle_displacements_, format_);
    char * data;
    particle_loader.read_msg(msg_initial, &data);
    particle_loader.copy_data_local(data);

  }

  if (sync_msg->next()) {
    // reset for next call (note not resetting at start since may get
    // called after messages received)
    sync_msg->reset();
    initialize_particle_mass(block);
    block->initial_done();
  }
  delete msg_initial;
}

//======================================================================

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
  cello::hierarchy()->root_blocks(root_blocks,  root_blocks+1, root_blocks+2);
  // Get this (reader) block index in root blocking
  index.array(array_lower, array_lower+1,array_lower+2);
  array_upper[0] = std::min(array_lower[0] + blocking_[0],root_blocks[0]);
  array_upper[1] = std::min(array_lower[1] + blocking_[1],root_blocks[1]);
  array_upper[2] = std::min(array_lower[2] + blocking_[2],root_blocks[2]);
}

void EnzoInitialHdf5::initialize_particle_mass(Block * block) {
  Particle particle = block->data()->particle();
  int divisor = std::pow(std::pow(2.0, cello::rank()), block->level());

  for (int it=0; it<particle.num_types(); it++) {
    if (particle.has_constant(it, "root_level_mass")) {
      //TODO: print warning if masses read from data file is being over written

      int ic = particle.constant_index(it, "root_level_mass");
      int ia = particle.attribute_index(it, "mass");

      union {
        float *  array_float;
        double * array_double;
      };

      int type_array = particle.attribute_type(it,ia);

      if (type_array == type_single) {
        float*  root_mass = (float*)  particle.constant_value(it, ic);
        initialize_particle_mass(array_float,  particle, it, ia, ((float)  *root_mass) / divisor);
      } else if (type_array == type_double) {
        double* root_mass = (double*) particle.constant_value(it, ic);
        initialize_particle_mass(array_double, particle, it, ia, ((double) *root_mass) / divisor);
      }
    }
  }
}

template <class T>
void EnzoInitialHdf5::initialize_particle_mass(T * array, Particle particle, int it, int ia, T mass) {
  int np = particle.num_particles(it);
  for (int ip=0; ip<np; ip++) {
    int ib, io;
    particle.index(ip,&ib,&io);
    array = (T*) particle.attribute_array(it,ia,ib);
    array[io] = mass;
  }
}






//=========================================================================
DataLoader::DataLoader(Block* block, std::string format) : block(block), m4()
{
  Field field = block->data()->field();
  // int index_field = field.field_id(name);
  // field.dimensions (0,&mx,&my,&mz); // (index_field,&mx,&my,&mz); int index_field = field.field_id(name);
  field.size         (&nx,&ny,&nz);

  block->lower(lower_block, lower_block+1, lower_block+2);
  block->upper(upper_block, upper_block+1, upper_block+2);
  format_ = format;
}

void DataLoader::open_file(std::string filename, std::string dataset, std::string coordinates) {
  file = new FileHdf5 ("./", filename);
  file->file_open();

  if (CHECK_COSMO_PARAMS)
    check_cosmology_(file);

  type_data = type_unknown;
  file->data_open(dataset, &type_data, m4, m4+1, m4+2, m4+3);

  ASSERT1("EnzoInitialHdf5::enforce_block()",
          "Unsupported type_data %d",
          type_data,
          ( (type_data == type_single) ||
            (type_data == type_double) ) );

  coords = coordinates;
}

void DataLoader::load(int* block_index, Index index_block) {
  char * data;
  read_dataset_(&data, index_block, block_index);

  (index_block == block->index()) ? copy_data_local(data) : copy_data_remote(index_block, data);
  delete_array_(&data, type_data);
}

void DataLoader::read_dataset_(char ** data, Index index_block, int block_index[3])
{
  // Get the grid size at level_
  // Hierarchy * hierarchy = cello::simulation()->hierarchy();
  // double lower_domain[3];
  // double upper_domain[3];
  // hierarchy->lower(lower_domain, lower_domain+1, lower_domain+2);
  // hierarchy->upper(upper_domain, upper_domain+1, upper_domain+2);

  // Read the domain dimensions
  IX = coords.find ("x");
  IY = coords.find ("y");
  IZ = coords.find ("z");

  ASSERT3 ("EnzoInitialHdf5::enforce_block()",
           "bad field coordinates %d %d %d",
           IX,IY,IZ,
           ((IX<4)&&(IY<4)&&(IZ<4)) &&
           ((IX != IY) || (IY==-1 && IZ == -1)) &&
           ((IX != IY && IY != IZ) || (IZ == -1)));

  // field size
  n4[0] = n4[1] = n4[2] = n4[3] = 1;
  n4[IX] = nx;
  n4[IY] = ny;
  n4[IZ] = nz;

  // compute cell widths
  h4[0] = h4[1] = h4[2] = h4[3] = 1.0;
  h4[IX] = (upper_block[0] - lower_block[0]) / (nx << index_block.level());
  h4[IY] = (upper_block[1] - lower_block[1]) / (ny << index_block.level());
  h4[IZ] = (upper_block[2] - lower_block[2]) / (nz << index_block.level());

  // determine offsets
  int o4[4] = {0,0,0,0};
  o4[IX] = block_index[0]*nx;
  o4[IY] = block_index[1]*ny;
  o4[IZ] = block_index[2]*nz;

  // open the dataspace
  file-> data_slice
    (m4[0],m4[1],m4[2],m4[3],
     n4[0],n4[1],n4[2],n4[3],
     o4[0],o4[1],o4[2],o4[3]);

  // create memory space
  // (fields was n4[(*IX)], n4[(*IY)], n4[(*IZ)])
  file->mem_create (nx,ny,nz,nx,ny,nz,0,0,0);

  // input domain size
  const int n = nx*ny*nz;
  (*data) = allocate_array_ (n, type_data);
  file->data_read ((*data));
}

void DataLoader::delete_array_(char ** array, int type_data)
{
  if (type_data == type_single) {
    delete [] (float*)(*array);
  } else if (type_data == type_double) {
    delete [] (double*)(*array);
  } else {
    ERROR1 ("DataLoader::delete_array_()",
            "Unsupported data type %d",type_data);
  }
  (*array) = nullptr;
}

char * DataLoader::allocate_array_ (int n, int type_data)
{
  char * data;
  if (type_data == type_single) {
    data = (char *)new float [n];
  } else if (type_data == type_double) {
    data = (char *)new double [n];
  } else {
    data = nullptr;
    ERROR1 ("DataLoader::allocate_array_()",
            "Unsupported data type %d",type_data);
  }
  return data;
}

void DataLoader::check_cosmology_(File * file) const
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




void FieldLoader::read_msg(MsgInitial * msg_initial, char ** data) {
  msg_initial->get_dataset(n4,h4,&nx,&ny,&nz,&IX,&IY,&IZ);
  msg_initial->get_field_data(&name, data, &type_data);
  Field field = block->data()->field();
  int index_field = field.field_id(name);
  field.dimensions (index_field,&mx,&my,&mz);
}

void FieldLoader::copy_data_local(char * data) {
  copy_dataset_to_field_(data);
}

void FieldLoader::copy_data_remote(Index index_block, char * data) {
  MsgInitial * msg_initial = new MsgInitial;
  msg_initial->set_dataset(n4,h4,nx,ny,nz,IX,IY,IZ);
  msg_initial->set_field_data(name, data, nx*ny*nz, type_data);
  enzo::block_array()[index_block].p_initial_hdf5_recv(msg_initial);
}

void FieldLoader::copy_dataset_to_field_(char * data) {
  Field field = block->data()->field();
  // Destination is this block--copy directly
  enzo_float * array = (enzo_float *) field.values(name);

  if (type_data == type_single) {
    copy_field_data_to_array_(array, (float *) data);
  } else if (type_data == type_double) {
    copy_field_data_to_array_(array, (double *) data);
  }
}

template <class T>
void FieldLoader::copy_field_data_to_array_(enzo_float * array, T * data) const {
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




void ParticleLoader::read_msg(MsgInitial * msg_initial, char ** data) {
  msg_initial->get_dataset (n4,h4,&nx,&ny,&nz,&IX,&IY,&IZ);
  int data_size;
  msg_initial->get_particle_data(&type,
                                 &attribute,
                                 data,&data_size,&type_data);
}

void ParticleLoader::copy_data_local(char * data) {
  copy_dataset_to_particle_(data);
}

void ParticleLoader::copy_data_remote(Index index_block, char * data) {
  MsgInitial * msg_initial = new MsgInitial;
  msg_initial->set_dataset (n4,h4,nx,ny,nz,IX,IY,IZ);
  msg_initial->set_particle_data(type, attribute,
                                data, nx*ny*nz, type_data);
  enzo::block_array()[index_block].p_initial_hdf5_recv(msg_initial);
}

void ParticleLoader::copy_dataset_to_particle_(char * data) {
  if (!block->is_leaf()) {return;}
  // Create particles and initialize them
  Particle particle = block->data()->particle();

  const int it = particle.type_index(type);
  const int ia = particle.attribute_index(it,attribute);
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
  block->lower(lower_block, lower_block+1, lower_block+2);

  if (l_particle_displacements_) {
    for (int axis=0; axis<3; axis++) {
      // if attribute corresponds to position along the given axis,
      // add to the attribute values the offsets of the cell centers
      // along that axis
      if (ia == particle.attribute_position(it,axis)) {
        if (type_array == type_single) {
          update_particle_displacements_
            (array_float,particle,it,ia,lower_block[axis],h4[IX],axis);
        } else {
          update_particle_displacements_
            (array_double,particle,it,ia,lower_block[axis],h4[IX],axis);
        }
      }
    }
  }
}

template <class T, class S>
void ParticleLoader::copy_particle_data_to_array_
(T * array, S * data, Particle particle, int it, int ia, int np)
{
  for (int ip=0; ip<np; ip++) {
    int ib,io;
    particle.index(ip,&ib,&io);
    array = (T*)particle.attribute_array(it,ia,ib);
    array[io] = data[ip];
  }
}

template <class T>
void ParticleLoader::update_particle_displacements_
(T * array, Particle particle, int it, int ia, double lower, double h, int axis)
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