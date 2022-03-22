// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodCheck.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2202-02-12
/// @brief    Implements the EnzoMethodCheck class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm.hpp"
#include "charm_enzo.hpp"
#include <iostream>
#include <iomanip>

// #define TRACE

#ifdef TRACE
#   undef TRACE
#   define TRACE(MSG)                                           \
  CkPrintf ("%d TRACE %s\n",CkMyPe(),std::string(MSG).c_str()); \
  fflush(stdout);
#   define TRACE_BLOCK(MSG,BLOCK)                               \
  CkPrintf ("%d TRACE %s %s\n",CkMyPe(),BLOCK->name().c_str(),  \
            std::string(MSG).c_str());                          \
  fflush(stdout);
#else
#   define TRACE(MSG)  /* ... */
#   define TRACE_BLOCK(MSG,BLOCK)  /* ... */
#endif

int Simulation::file_counter_ = 0;

//----------------------------------------------------------------------

EnzoMethodCheck::EnzoMethodCheck
(int num_files, std::string ordering, std::vector<std::string> directory)
  : Method(),
    num_files_(num_files),
    ordering_(ordering),
    directory_(directory)
{
  TRACE("[1] EnzoMethodCheck::EnzoMethodCheck()");
  Refresh * refresh = cello::refresh(ir_post_);
  cello::simulation()->refresh_set_name(ir_post_,name());
  refresh->add_field("density");
  // Create IO writer
  if (CkMyPe() == 0) {

    enzo::simulation()->set_sync_check_writer(num_files_);

    proxy_io_enzo_writer = CProxy_IoEnzoWriter::ckNew
      (num_files, ordering, num_files);

    proxy_io_enzo_writer.doneInserting();

    proxy_enzo_simulation.p_set_io_writer(proxy_io_enzo_writer);

  }

}

//----------------------------------------------------------------------

void EnzoSimulation::p_set_io_writer(CProxy_IoEnzoWriter io_enzo_writer)
{ proxy_io_enzo_writer = io_enzo_writer; }

//----------------------------------------------------------------------

void EnzoMethodCheck::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | num_files_;
  p | ordering_;
  p | directory_;

}

//----------------------------------------------------------------------

void EnzoMethodCheck::compute ( Block * block) throw()
{
  TRACE_BLOCK("[2] EnzoMethodCheck::compute()",block);
  CkCallback callback(CkIndex_EnzoSimulation::r_method_check_enter(NULL),0,
                      proxy_enzo_simulation);
  block->contribute(callback);
}

//----------------------------------------------------------------------

void EnzoSimulation::r_method_check_enter(CkReductionMsg *msg)
// [ Called on ip=0 only ]
{
  TRACE("[3] EnzoSimulation::r_method_check_enter()");
  
  delete msg;

  check_num_files_  = enzo::config()->method_check_num_files;
  check_ordering_   = enzo::config()->method_check_ordering;
  check_directory_  = enzo::config()->method_check_dir;

  /// Initialize synchronization counters
  sync_check_done_.          set_stop(check_num_files_);

  /// Create the directory

  bool already_exists = false;
  
  std::string name_dir = file_create_dir_(check_directory_,already_exists);

  if (already_exists) {
    // Exit checkpoint if directory already exists
    enzo::block_array().p_check_done();
    
  } else {
    // Else start checkpoint

    // Create hierarchy file if root writer

    std::string name_file = name_dir + "/check.file_list";
    std::ofstream stream_file_list (name_file);

    ASSERT1("IoEnzoWriter",
            "Cannot open hierarchy file %s for writing",
            name_file,stream_file_list);

    stream_file_list << std::setfill('0');
    int max_digits = log(check_num_files_-1)/log(10) + 1;
    stream_file_list << check_num_files_ << "\n";
    for (int i=0; i<check_num_files_; i++) {
      stream_file_list << "block_data-" << std::setw(max_digits) << i << "\n";
    }
    stream_file_list.flush();

    enzo::block_array().p_check_write_first
      (check_num_files_, check_ordering_, name_dir);
  }
  // Create IoEnzoWriter array. Synchronizes by calling
  // EnzoSimulation[0]::p_writer_created() when done

}

//----------------------------------------------------------------------

IoEnzoWriter::IoEnzoWriter
(int num_files, std::string ordering) throw ()
  : CBase_IoEnzoWriter(),
    num_files_(num_files),
    ordering_(ordering)
{
  TRACE("[4] IoEnzoWriter::IoEnzoWriter()");
}

//----------------------------------------------------------------------

void EnzoBlock::p_check_write_first
(int num_files, std::string ordering, std::string name_dir)
{
  TRACE_BLOCK("[8] EnzoBlock::p_check_write_first",this);

  EnzoMsgCheck * msg_check;
  bool is_first (false);
  const int index_file = create_msg_check_
    (&msg_check,num_files,ordering,name_dir,&is_first);

  if (is_first) {
    proxy_io_enzo_writer[index_file].p_write (msg_check);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_check_write_next(int num_files, std::string ordering)
{
  TRACE_BLOCK("[9] EnzoBlock::p_check_write_next",this);

  std::string name_dir {""};
  EnzoMsgCheck * msg_check;
  const int index_file = create_msg_check_
    (&msg_check,num_files,ordering);

  proxy_io_enzo_writer[index_file].p_write (msg_check);
}

//----------------------------------------------------------------------

void IoEnzoWriter::p_write (EnzoMsgCheck * msg_check)
{
  std::string name_this, name_next;
  Index index_this, index_next;
  int index_block;
  bool is_first, is_last;
  std::string name_dir;

  msg_check->get_parameters
    (index_this,index_next,name_this,name_next,
     index_block,is_first,is_last,name_dir);

  // Write to block list file, opening or closing file as needed

  CkPrintf ("TRACE_CHECK recv p_write first %d last %d name %s\n",
            is_first?1:0,is_last?1:0,name_this.c_str());
  
  if (is_first) {
    // Create block list
    stream_block_list_ = create_block_list_(name_dir);

    // Create HDF5 file
    char name_file[80];
    sprintf (name_file,"block_data-%d.h5", thisIndex);
    file_ = file_open_(name_dir,name_file);

    // Write HDF5 header meta data
    file_write_hierarchy_();
  }

  // Write block list
  write_block_list_(name_this);

  // Write Block to HDF5
  file_write_block_(msg_check);

  if (is_last) {
    // close block list
    close_block_list_();
    // close HDF5 file
    file_->file_close();
  }

  TRACE("[A] IoEnzoWriter::p_write_first");
  if (!is_last) {
    enzo::block_array()[index_next].p_check_write_next(num_files_, ordering_);
  } else {
    proxy_enzo_simulation[0].p_check_done();
  }
}

//----------------------------------------------------------------------

std::ofstream IoEnzoWriter::create_block_list_(std::string name_dir)
{
  char file_name[80];
  sprintf (file_name,"%s/block_data-%d.block_list",
           name_dir.c_str(),thisIndex);

  std::ofstream stream_block_list (file_name);

  ASSERT1("Simulation::create_block_list_",
          "Cannot open block_list file %s for writing",
          file_name,stream_block_list);

  return stream_block_list;
}

//----------------------------------------------------------------------

void IoEnzoWriter::write_block_list_(std::string block_name)
{
  stream_block_list_ << block_name << "\n";
}

//----------------------------------------------------------------------

void IoEnzoWriter::close_block_list_()
{
  stream_block_list_.flush();
}

//----------------------------------------------------------------------

void EnzoSimulation::p_check_done()
{
  TRACE("[B] EnzoSimulation::p_check_done");
  if (sync_check_done_.next()) {
    enzo::block_array().p_check_done();
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_check_done()
{
  TRACE_BLOCK("[C] EnzoBlock::p_check_done",this);
  compute_done();
}

//======================================================================

int EnzoBlock::create_msg_check_
( EnzoMsgCheck ** msg_check,
  int num_files, std::string ordering,
  std::string name_dir,
  bool *          is_first
  )
{
  ScalarData<int> *   scalar_data_int    = data()->scalar_data_int();
  ScalarDescr *       scalar_descr_int   = cello::scalar_descr_int();
  ScalarData<Index> * scalar_data_index  = data()->scalar_data_index();
  ScalarDescr *       scalar_descr_index = cello::scalar_descr_index();
  const int is_index = scalar_descr_int->index(ordering+":index");
  const int is_count = scalar_descr_int->index(ordering+":count");
  const int is_next  = scalar_descr_index->index(ordering+":next");
  const int count    = *scalar_data_int->value(scalar_descr_int,is_count);
  const Index next   = *scalar_data_index->value(scalar_descr_index,is_next);

  const int rank = cello::rank();
  const Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  int na3[3];
  hierarchy->root_blocks(na3,na3+1,na3+2);

  const int min_level = hierarchy->min_level();

  Index index_this, index_next;
  std::string name_this, name_next;
  int index_block, index_file;
  bool is_last;

  index_this = this->index();
  index_next = next;

  name_this = this->name();
  name_next = this->name(index_next);

  index_block    = *scalar_data_int->value(scalar_descr_int,is_index);
  index_file = index_block*num_files/count;

  const int ib  = index_block;
  const int ibm = index_block - 1;
  const int ibp = index_block + 1;
  const int nb = count;
  const int nf = num_files;
  if (is_first) { (*is_first) = (ib  == 0)  || (ib*nf/nb != (ibm)*nf/nb); }
  is_last  = (ibp == nb) || (ib*nf/nb != (ibp)*nf/nb);

  *msg_check = new EnzoMsgCheck;

  (*msg_check)->set_block(this);

  (*msg_check)->set_parameters
    (index_this,index_next,name_this,name_next,
     index_block,is_first?(*is_first):false,is_last);

  (*msg_check)->set_name_dir (name_dir);

  return index_file;
}

//----------------------------------------------------------------------

std::string Simulation::file_create_dir_
(std::vector<std::string> directory_format, bool & already_exists)
{
  const int counter = Simulation::file_counter_++;
  const int cycle = cello::simulation()->cycle();
  const double time = cello::simulation()->time();

  cello::create_directory
    (&directory_format, counter,cycle,time,already_exists);

  ASSERT("Simulation::file_create_dir_",
         "Directory name must be non-empty",
         directory_format[0] != "");

  std::string name_dir = cello::expand_name
    (&directory_format,counter, cycle, time);

  return name_dir;
}

//----------------------------------------------------------------------
FileHdf5 * IoEnzoWriter::file_open_
(std::string path_name, std::string file_name)
{
  // Create File
  FileHdf5 * file = new FileHdf5 (path_name, file_name);
  file->file_create();

  return file;
}

//----------------------------------------------------------------------

void IoEnzoWriter::file_write_hierarchy_()
{
  IoHierarchy io_hierarchy = (cello::hierarchy());
  for (size_t i=0; i<io_hierarchy.meta_count(); i++) {

    void * buffer;
    std::string name;
    int type_scalar;
    int nx,ny,nz;

    // Get object's ith metadata
    io_hierarchy.meta_value(i,& buffer, &name, &type_scalar, &nx,&ny,&nz);

    // Write object's ith metadata
    file_->file_write_meta(buffer,name.c_str(),type_scalar,nx,ny,nz);
  }
}

//----------------------------------------------------------------------

void IoEnzoWriter::file_write_block_ (EnzoMsgCheck * msg_check)
{
  Index  index_block;
  Index  index_next;
  std::string  name_block;
  std::string  name_next;
  int  block_order;
  bool  is_first;
  bool  is_last;
  std::string  name_dir;

  msg_check->get_parameters
    ( index_block,
      index_next,
      name_block,
      name_next,
      block_order,
      is_first,
      is_last,
      name_dir);

  //  const bool is_local = (msg_check == nullptr);

  IoBlock * io_block = msg_check->io_block();
  CkPrintf ("TRACE_CHECK io_block = %p\n",(void*)io_block);
  fflush(stdout);
  double * lower = msg_check->block_lower();
  double * upper = msg_check->block_upper();
  int * size     = msg_check->block_size();

  // Create file group for block

  std::string group_name = "/" + name_block;

  file_->group_chdir(group_name);
  file_->group_create();

  // Write block meta data

  write_meta_ (file_, io_block, "group");

  // // Create new data object to hold EnzoMsgCheck/DataMsg fields and particles

  Data * data;
  bool data_allocated (true);

  // @@@@ DataMsg does not currently support num_field_data > 1
  
  int num_field_data=1;
  data = new Data
    (size[0],size[1],size[2],
     num_field_data,
     lower[0],lower[1],lower[2],
     upper[0],upper[1],upper[2],
     cello::field_descr(),
     cello::particle_descr());

  data->allocate();

  msg_check->update(data);

  // Write Block Field data

  // number of "history" field data objects
  const int nh = data->num_field_data();
  const int nf = cello::field_descr()->field_count();
  // // May have multiple field_data objects for field history
  if (nh > 0) {
    for (int i_h=0; i_h<nh; i_h++) {
      FieldData * field_data = data->field_data(i_h);
      for (int i_f=0; i_f<nf; i_f++) {
        const int index_field = i_f;
        IoFieldData * io_field_data = enzo::factory()->create_io_field_data();

        void * buffer;
        std::string name;
        int type;
        int mx,my,mz;  // Array dimension
        int nx,ny,nz;  // Array size

        io_field_data->set_field_data((FieldData*)field_data);
        io_field_data->set_field_index(index_field);

        io_field_data->field_array
          (&buffer, &name, &type, &mx,&my,&mz, &nx,&ny,&nz);

        file_->mem_create(nx,ny,nz,nx,ny,nz,0,0,0);
        if (mz > 1) {
          file_->data_create(name.c_str(),type,mz,my,mx,1,nz,ny,nx,1);
        } else if (my > 1) {
          file_->data_create(name.c_str(),type,my,mx,  1,1,ny,nx, 1,1);
        } else {
          file_->data_create(name.c_str(),type,mx,  1,  1,1,nx,  1,1,1);
        }
        file_->data_write(buffer);
        file_->data_close();

        delete io_field_data;
      }
    }
  }
  // // Write Block Particle data

  Particle particle = data->particle();

  for (int it=0; it<particle.num_types(); it++) {

    // get the number of particle batches and attributes
    const int nb = particle.num_batches(it);
    const int na = particle.num_attributes(it);

    // For each particle attribute
    for (int ia=0; ia<na; ia++) {

      IoParticleData * io_particle_data =
        enzo::factory()->create_io_particle_data();

      // For each particle attribute
      int np = particle.num_particles (it);

      const std::string name = "particle_"
        +                particle.type_name(it) + "_"
        +                particle.attribute_name(it,ia);

      const int type = particle.attribute_type(it,ia);

      // create the disk array
      file_->data_create(name.c_str(),type,np,1,1,1,np,1,1,1);

      // running count of particles in the type
      int i0 = 0;

      // for each batch of particles
      for (int ib=0; ib<nb; ib++) {

        // number of particles in batch (it,ib)
        const int mb = particle.num_particles(it,ib);

        // create the memory space for the batch
        file_->mem_create(mb,1,1,mb,1,1,0,0,0);

        // get the buffer of data
        const void * buffer = (const void *) particle.attribute_array(it,ia,ib);

        // find the hyper_slab of the disk dataset
        file_->data_slice (np,1,1,1, mb,1,1,1, i0,0,0,0);

        // update the running count of particles for the type
        i0 += mb;

        // write the batch to disk
        file_->data_write(buffer);

        // close the memory space
        file_->mem_close();
      }

      // check that the number of particles equals the number written

      ASSERT2 ("OutputData::write_particle_data()",
               "Particle count mismatch %d particles %d written",
               np,i0,
               np == i0);

      // close the attribute dataset
      file_->data_close();
      delete io_particle_data;
    }
  }

  if (data_allocated) delete data;

  file_->group_close();
}

//----------------------------------------------------------------------

DataMsg * EnzoMethodCheck::create_data_msg_ (Block * block)
{
  int if3[3] = {0,0,0};
  int ic3[3] = {0,0,0};
  int g3[3] = {0};

  // Create refresh object required by FieldFace
  Refresh * refresh = new Refresh;

  // Initialize refresh fields
  bool any_fields = false;
  if (cello::field_descr()->field_count() > 0) {
    refresh->add_all_fields();
    any_fields = true;
  }

  // Initialize refresh particles
  bool any_particles = false;
  if (cello::particle_descr()->num_types()) {
    refresh->add_all_particles();
    any_particles = true;
  }

  // Create FieldFace object specifying fields to send
  FieldFace * field_face = block->create_face
    (if3,ic3,g3, refresh_same, refresh, true);

  // Create data message object to send
  DataMsg * data_msg = new DataMsg;
  if (any_fields) {
    data_msg -> set_field_face (field_face,true);
    data_msg -> set_field_data (block->data()->field_data(),false);
  }
  if (any_particles) {
    data_msg -> set_particle_data (block->data()->particle_data(),false);
  }
  return data_msg;
}

//----------------------------------------------------------------------

void IoEnzoWriter::write_meta_
( FileHdf5 * file, Io * io, std::string type_meta )
{
  for (size_t i=0; i<io->meta_count(); i++) {

    void * buffer;
    std::string name;
    int type_scalar;
    int nx,ny,nz;

    // Get object's ith metadata

    io->meta_value(i,& buffer, &name, &type_scalar, &nx,&ny,&nz);

    // Write object's ith metadata
    if ( type_meta == "group" ) {
      file->group_write_meta(buffer,name.c_str(),type_scalar,nx,ny,nz);
    } else if (type_meta == "file") {
      file->file_write_meta(buffer,name.c_str(),type_scalar,nx,ny,nz);
    } else {
      ERROR1 ("MethodOutput::write_meta_()",
              "Unknown type_meta \"%s\"",
              type_meta.c_str());
    }
  }
}

