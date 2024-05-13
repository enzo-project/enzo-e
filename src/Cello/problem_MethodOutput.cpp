// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOutput.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-04-27
/// @brief    Implementation of the Output "method"

#include "problem.hpp"
#include "charm_simulation.hpp"
#include "test.hpp"

// #define TRACE_OUTPUT

// #define BYPASS_BLOCK_TRACE

//----------------------------------------------------------------------

// this is a commonly occuring operation that should probably be directly
// supported by ParameterGroup
static std::vector<std::string> read_str_vec_(ParameterGroup p,
                                              std::string name) noexcept
{
  int length = p.list_length(name);
  std::vector<std::string> out(length);
  for (int i = 0; i < length; i++) { out[i] = p.list_value_string(i, name); }
  return out;
}

//----------------------------------------------------------------------

static std::vector<std::string> read_path_param_(ParameterGroup p,
                                                 std::string name) noexcept
{
  const std::string root_path = p.get_group_path();
  if (p.type(name) == parameter_string) {
    return { p.value_string(name,"") };
  } else if (p.type(name) == parameter_list) {
    std::vector<std::string> out = read_str_vec_(p, name);
    if (out.size() == 0){
      ERROR2("read_path_param_", "%s:%s is empty",
             root_path.c_str(), name.c_str());
    }
    return out;
  }
  ERROR2("read_path_param_", "%s:%s must be a string or a list of strings",
         root_path.c_str(), name.c_str());
}

//----------------------------------------------------------------------

MethodOutput::MethodOutput(const Factory * factory, ParameterGroup p) noexcept
  : MethodOutput(factory,
                 read_path_param_(p, "file_name"),
                 read_path_param_(p, "path_name"),
                 read_str_vec_(p, "field_list"),
                 read_str_vec_(p, "particle_list"),
                 p.value_integer("ghost_depth",0),
                 p.value_integer("min_face_rank",0), // default 0 all faces
                 p.value_logical("all_fields", false),
                 p.value_logical("all_particles", false),
                 p.value_logical("all_blocks", true),
                 p.list_value_integer(0,"blocking",1),
                 p.list_value_integer(1,"blocking",1),
                 p.list_value_integer(2,"blocking",1))
{ }

//----------------------------------------------------------------------

MethodOutput::MethodOutput
  (const Factory * factory,
   std::vector< std::string > file_name,
   std::vector< std::string > path_name,
   std::vector< std::string > field_list,
   std::vector< std::string > particle_list,
   int ghost_depth,
   int min_face_rank,
   bool all_fields,
   bool all_particles,
   bool all_blocks,
   int blocking_x,
   int blocking_y,
   int blocking_z) noexcept
    : Method(),
      file_name_(file_name),
      path_name_(path_name),
      field_list_(),
      particle_list_(),
      ghost_depth_(ghost_depth),
      min_face_rank_(min_face_rank),
      all_fields_(all_fields),
      all_particles_(all_particles),
      blocking_(),
      is_count_(-1),
      is_block_list_(-1),
      factory_(factory),
      all_blocks_(all_blocks)
{
  if (field_list.size() > 0) {
    field_list_.resize(field_list.size());
    for (size_t i=0; i<field_list.size(); i++) {
      const int index_field = cello::field_descr()->field_id(field_list[i]);
      field_list_[i] = index_field;
    }
  }

  if (particle_list.size() > 0) {
    particle_list_.resize(particle_list.size());
    for (size_t i=0; i<particle_list.size(); i++) {
      const int index_particle =
        cello::particle_descr()->type_index(particle_list[i]);
      particle_list_[i] = index_particle;
    }
  }

  Refresh * refresh = cello::refresh(ir_post_);

  refresh->set_ghost_depth (ghost_depth);
  refresh->set_min_face_rank(min_face_rank);

  // add fields to refresh
  if (all_fields_) {
    refresh->add_all_fields();
  } else {
    const int nf=field_list_.size();
    for (int i_f=0; i_f<nf; i_f++) {
      cello::refresh(ir_post_)->add_field(field_list_[i_f]);
    }
  }

  // add particles to refresh
  if (all_particles_) {
    refresh->add_all_particles();
  } else {
    const int nf=particle_list_.size();
    for (int i_f=0; i_f<nf; i_f++) {
      cello::refresh(ir_post_)->add_particle(particle_list_[i_f]);
    }
  }

  // Create field list if all_fields_
  if (all_fields_ && field_list_.size() == 0) {
    int nf = cello::field_descr()->field_count();
    for (int i_f=0; i_f<nf; i_f++) {
      field_list_.push_back(i_f);
    }
  }

  // Create particle list if all_particles_
  if (all_particles_ && particle_list_.size() == 0) {
    int np = cello::particle_descr()->num_types();
    for (int i_p=0; i_p<np; i_p++) {
      particle_list_.push_back(i_p);
    }
  }

  ASSERT("MethodOutput()",
         "MethodOutput requires either field_list or particle_list "
         "to be non-empty",
         (all_fields || all_particles ||
          field_list.size() > 0 || particle_list.size() > 0));

  // initialize blocking partition
  blocking_[0] = blocking_x;
  blocking_[1] = blocking_y;
  blocking_[2] = blocking_z;

  ScalarDescr * sdi = cello::scalar_descr_int();
  is_count_ = sdi->new_value("method_output:count");

  ScalarDescr * sdp = cello::scalar_descr_void();
  is_block_list_ = sdp->new_value("method_output:block_list");
}

//----------------------------------------------------------------------

MethodOutput::~MethodOutput() throw()
{
}

//----------------------------------------------------------------------

void MethodOutput::pup (PUP::er &p)
{
  TRACEPUP;
  Method::pup(p);
  p | file_name_;
  p | path_name_;
  p | field_list_;
  p | particle_list_;
  p | ghost_depth_;
  p | min_face_rank_;
  p | all_fields_;
  p | all_particles_;
  PUParray(p,blocking_,3);
  p | is_count_;
  p | is_block_list_;
  p | factory_nonconst_;
  p | all_blocks_;
}

//----------------------------------------------------------------------

void MethodOutput::compute ( Block * block) throw()
{
  // barrier to ensure tree traversal doesn't reach a block
  // before the method has started

  if (block->index().is_root()) {
    // Create directory if any
    ScalarData<int> * scalar_int = block->data()->scalar_data_int();
    int * counter = scalar_int->value(cello::scalar_descr_int(),is_count_);
    (*counter)++;
    bool already_exists = false;
    std::string path_name = cello::create_directory
      (&path_name_,(*counter),block->cycle(),block->time(),already_exists);
  }

  CkCallback callback(CkIndex_Block::r_method_output_continue(nullptr),
                      cello::block_array());
  block->contribute(callback);

}

//----------------------------------------------------------------------

void Block::r_method_output_continue(CkReductionMsg *msg)
{
  delete msg;
  MethodOutput * method = static_cast<MethodOutput*> (this->method());
  method->compute_continue(this);
}

//----------------------------------------------------------------------

void MethodOutput::compute_continue(Block * block)
{
  // called by Block::r_method_output_continue()

  // negative level blocks do not participate...
  if (block->level() < 0) {
    compute_done(block);
    return;
  }
  // non-writers wait until called later...
  if (! is_writer_(block->index()) ) {
    return;
  }

  // otherwise this block is a writer, so start writing
  int a3[3];
  block->index().array(a3,a3+1,a3+2);

  int root_blocks[3];
  cello::hierarchy()->root_blocks
    (root_blocks,root_blocks+1,root_blocks+2);

  int index_min[3] = {a3[0],a3[1],a3[2]};
  int index_max[3] = {std::min(a3[0] + blocking_[0],root_blocks[0]),
                      std::min(a3[1] + blocking_[1],root_blocks[1]),
                      std::min(a3[2] + blocking_[2],root_blocks[2])};

  BlockTrace bt (cello::rank(),index_min,index_max);

  FileHdf5 * file = file_open_(block,a3);

  // Open *.block_list file and save FILE pointer
  const int count   = file_count_(block);
  const int cycle   = block->cycle();
  const double time = block->time();
  std::string file_name = cello::expand_name
    (&file_name_,count,cycle,time);

  ScalarData<void *> * scalar_void = block->data()->scalar_data_void();
  FILE ** fp_block_list =
    (FILE **) scalar_void->value(cello::scalar_descr_void(),is_block_list_);
  *fp_block_list = fopen ((file_name+".block_list").c_str(),"w");

  // Write and close file_list file
  FILE * fp_file_list = fopen ((file_name+".file_list").c_str(),"w");
  fprintf (fp_file_list,"%s\n",file_name.c_str());
  fclose(fp_file_list);

  // Restore working directory now that files are created
  // (changed for block_list and file_list files in file_open_())
  if (chdir ("..") != 0) {
    ERROR("MethodOutput::compute_continue",
          "error occured while changing current directory to parent dir");
  }

  // Create output message
  MsgOutput * msg_output = new MsgOutput(bt,this,file);

  msg_output->set_index_send (block->index());
  // Get its block_trace object
  BlockTrace * block_trace = msg_output->block_trace();

  // write the version number to file
  cello::io::write_version_metadata(file);

  // write hierarchy meta-data to file
  file_write_hierarchy_(file);

  if (all_blocks_ || block->is_leaf()) {
    // write this (writer) block's data to file
    msg_output->set_block(block,factory_);
    file_write_block_(file,block,nullptr);

  }

  // Check if any non-writer blocks
  if ( ! block_trace->next(block->is_leaf())) {
    Index index_next = block_trace->top();
    msg_output->del_block();
    // Initiate traversing distributed octree forest to request
    // blocks to send their data to this writer block to write
    cello::block_array()[index_next].p_method_output_next(msg_output);

  } else {
    // Close file
    FileHdf5 * file = msg_output->file();
    file->file_close();
    delete file;

    // Close *.block_list file
    ScalarData<void *> * scalar_void = block->data()->scalar_data_void();
    FILE ** fp_block_list = (FILE **)
      scalar_void->value(cello::scalar_descr_void(),is_block_list_);
    fclose(*fp_block_list);
    *fp_block_list = nullptr;

    // Delete message
    delete msg_output;
    msg_output = nullptr;

    // Exit MethodOutput
    compute_done(block);
  }
}

//----------------------------------------------------------------------

void Block::p_method_output_next (MsgOutput * msg_output)
{
  MethodOutput * method_output =
    static_cast<MethodOutput*>(cello::problem()->method(index_method_));
  method_output -> next(this,msg_output);
}

//----------------------------------------------------------------------

void Block::p_method_output_write (MsgOutput * msg_output)
{
  MethodOutput * method_output =
    static_cast<MethodOutput*>(cello::problem()->method(index_method_));
  method_output -> write(this,msg_output);
}

//----------------------------------------------------------------------

void MethodOutput::next(Block * block, MsgOutput * msg_output_in )
// Process the next block in the depth-first distributed octree forest
// traversal.  Send data to writer if a leaf, or go directly to the
// next Block in the traversal if not
{

  // copy incoming message and delete old (cannot reuse!)
  MsgOutput * msg_output = new MsgOutput (*msg_output_in);

  delete msg_output_in;
  msg_output_in = nullptr;

  const bool is_leaf = block->is_leaf();
  BlockTrace * bt = msg_output->block_trace();
  bt->next(block->is_leaf());
  msg_output->set_index_send (block->index());

  if (all_blocks_ || is_leaf) {
    // if leaf (or writing all blocks), copy data to msg and send to writer
#ifdef BYPASS_BLOCK_TRACE
    Index index_home = block->index().next(rank,na3,block->is_leaf());
#else
    Index index_home = bt->home();
#endif
    DataMsg * data_msg = create_data_msg_(block);
    msg_output->set_data_msg(data_msg);
    msg_output->set_block(block,factory_);
    cello::block_array()[index_home].p_method_output_write(msg_output);
  } else {
    // if non-leaf (and not all blocks), forward message to child
    Index index_next = bt->top();
    msg_output->del_block();
    cello::block_array()[index_next].p_method_output_next(msg_output);
  }
  compute_done(block);
}

//----------------------------------------------------------------------

void MethodOutput::write(Block * block, MsgOutput * msg_output_in )
// Writes incoming data from the distributed octree forest traversal,
// and continues to the next block in the distributed octree forest
// traversal (if any)
{
  // Write the block data
  FileHdf5 * file = msg_output_in->file();
  file_write_block_(file,block,msg_output_in);

  // Copy incoming message
  MsgOutput * msg_output = new MsgOutput (*msg_output_in);
  // and delete it (cannot reuse)
  delete msg_output_in;
  msg_output_in = nullptr;

  msg_output->set_index_send (block->index());
  BlockTrace * bt = msg_output->block_trace();
  Index index_next = bt->top();
  Index index_home = bt->home();
  if (index_next == index_home) {
    // done
    // Close file
    FileHdf5 * file = msg_output->file();
    file->file_close();
    delete file;
    delete msg_output;

    // Close *.block_list file
    ScalarData<void *> * scalar_void = block->data()->scalar_data_void();
    FILE ** fp_block_list = (FILE **)
      scalar_void->value(cello::scalar_descr_void(),is_block_list_);
    fclose(*fp_block_list);
    *fp_block_list = nullptr;

    compute_done(block);
  } else {
    msg_output->del_block();
    cello::block_array()[index_next].p_method_output_next(msg_output);
  }
}

//----------------------------------------------------------------------

void MethodOutput::compute_done (Block * block)
{
  CkCallback callback(CkIndex_Block::r_method_output_done(nullptr), 
                      cello::block_array());
  block->contribute(callback);
}

//----------------------------------------------------------------------

void Block::r_method_output_done(CkReductionMsg *msg)
{
  delete msg;
  this->compute_done();
}

//======================================================================

int MethodOutput::is_writer_ (Index index)
{
  int a3[3];
  index.array(a3,a3+1,a3+2);
  int level = index.level();
  return ((level==0) &&
          ( a3[0] % blocking_[0] == 0) &&
          ( a3[1] % blocking_[1] == 0) &&
          ( a3[2] % blocking_[2] == 0));
}

//----------------------------------------------------------------------
FileHdf5 * MethodOutput::file_open_(Block * block, int a3[3])
{
  // Generate file path
  ScalarData<int> * scalar_int = block->data()->scalar_data_int();
  int * counter = scalar_int->value(cello::scalar_descr_int(),is_count_);

  std::string path_name = cello::expand_name
    (&path_name_,(*counter),block->cycle(),block->time());

  // Generate file name
  int count = file_count_(block);
  std::string file_name = cello::expand_name
    (&file_name_,count,block->cycle(),block->time());

  if (block->index().is_root()) {
    Monitor::instance()->print
      ("Output","MethodOutput writing data file %s",
       (path_name + "/" + file_name).c_str());
  }

  // Create File
  FileHdf5 * file = new FileHdf5 (path_name, file_name);
  file->file_create();

  // Change directory for file_list and block_list files
  if (chdir (path_name.c_str()) != 0) {
    ERROR1("MethodOutput::file_open_",
           "Error occured while changing directory to \"%s\"",
           path_name.c_str());
  }

  return file;
}

//----------------------------------------------------------------------

void MethodOutput::file_write_hierarchy_(FileHdf5 * file)
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
    file->file_write_meta(buffer,name.c_str(),type_scalar,nx,ny,nz);
  }
}

//----------------------------------------------------------------------

void MethodOutput::file_write_block_
(FileHdf5 * file, Block * block, MsgOutput * msg_output)
{
  const bool is_local = (msg_output == nullptr);

  std::string block_name = (is_local) ?
    block->name() : msg_output->block_name();

  IoBlock * io_block = (is_local) ?
    factory_->create_io_block() : msg_output->io_block();
  if (is_local) io_block->set_block(block);

  // Write Block name to block_list file
  ScalarData<void *> * scalar_void = block->data()->scalar_data_void();
  FILE ** fp_block_list = (FILE **)
    scalar_void->value(cello::scalar_descr_void(),is_block_list_);
  fprintf(*fp_block_list,"%s\n",block_name.c_str());

  double block_lower[3];
  double block_upper[3];
  if (is_local) {
    block->lower(block_lower,block_lower+1,block_lower+2);
    block->upper(block_upper,block_upper+1,block_upper+2);
  }

  double * lower = (is_local)?
    block_lower : msg_output->block_lower();
  double * upper = (is_local)?
    block_upper : msg_output->block_upper();

  // Create file group for block

  std::string group_name = "/" + block_name;

  file->group_chdir(group_name);
  file->group_create();

  // Write block meta data

  write_meta_ (file, io_block, "group");

  // Create new data object to hold MsgOutput/DataMsg fields and particles

  int nx,ny,nz;
  int num_field_data = 1;
  block->data()->field().size(&nx,&ny,&nz);

  Data * data;

  bool data_allocated;
  if (is_local) {
    data_allocated = false;
    data = block->data();
  } else {
    data_allocated = true;
    data = new Data
    (nx,ny,nz,
     num_field_data,
     lower[0],lower[1],lower[2],
     upper[0],upper[1],upper[2],
     cello::field_descr(),
     cello::particle_descr());

    // Allocate fields
    data->allocate();

    msg_output->update(data);
  }

  // Write Block Field data

  FieldData * field_data = data->field_data();

  for (std::size_t i_f=0; i_f<field_list_.size(); i_f++) {
    const int index_field = field_list_[i_f];
    IoFieldData * io_field_data = factory_->create_io_field_data();

    void * buffer;
    std::string name;
    int type;
    int mx,my,mz;  // Array dimension
    int nx,ny,nz;  // Array size

    io_field_data->set_field_data((FieldData*)field_data);
    io_field_data->set_field_index(index_field);

    io_field_data->field_array
      (&buffer, &name, &type, &mx,&my,&mz, &nx,&ny,&nz);

    file->mem_create(nx,ny,nz,nx,ny,nz,0,0,0);
    if (mz > 1) {
      file->data_create(name.c_str(),type,mz,my,mx,1,nz,ny,nx,1);
    } else if (my > 1) {
      file->data_create(name.c_str(),type,my,mx,  1,1,ny,nx, 1,1);
    } else {
      file->data_create(name.c_str(),type,mx,  1,  1,1,nx,  1,1,1);
    }
    file->data_write(buffer);
    file->data_close();

    delete io_field_data;
  }

  // Write Block Particle data

  Particle particle = data->particle();

  for (std::size_t i_p=0; i_p<particle_list_.size(); i_p++) {

    // Get particle type for it'th element of the particle output list
    const int it = particle_list_[i_p];

    // get the number of particle batches and attributes
    const int nb = particle.num_batches(it);
    const int na = particle.num_attributes(it);

    // For each particle attribute
    for (int ia=0; ia<na; ia++) {

      IoParticleData * io_particle_data = factory_->create_io_particle_data();

      // For each particle attribute
      int np = particle.num_particles (it);

      const std::string name = "particle_"
        +                particle.type_name(it) + "_"
        +                particle.attribute_name(it,ia);

      const int type = particle.attribute_type(it,ia);

      // create the disk array
      file->data_create(name.c_str(),type,np,1,1,1,np,1,1,1);

      // running count of particles in the type
      int i0 = 0;

      // for each batch of particles
      for (int ib=0; ib<nb; ib++) {

        // number of particles in batch (it,ib)
        const int mb = particle.num_particles(it,ib);

        // create the memory space for the batch
        file->mem_create(mb,1,1,mb,1,1,0,0,0);

        // get the buffer of data
        const void * buffer = (const void *) particle.attribute_array(it,ia,ib);

        // find the hyper_slab of the disk dataset
        file->data_slice (np,1,1,1, mb,1,1,1, i0,0,0,0);

        // update the running count of particles for the type
        i0 += mb;

        // write the batch to disk
        file->data_write(buffer);

        // close the memory space
        file->mem_close();
      }

      // check that the number of particles equals the number written

      ASSERT2 ("OutputData::write_particle_data()",
               "Particle count mismatch %d particles %d written",
               np,i0,
               np == i0);

      // close the attribute dataset
      file->data_close();
      delete io_particle_data;
    }
  }

  if (data_allocated) delete data;
  file->group_close();
}

//----------------------------------------------------------------------

void MethodOutput::write_meta_ ( FileHdf5 * file, Io * io, std::string type_meta )
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

//----------------------------------------------------------------------

int MethodOutput::file_count_(Block * block)
{
  int ax,ay,az;
  cello::hierarchy()->root_blocks(&ax,&ay,&az);
  int mx=(ax+blocking_[0]-1)/blocking_[0];
  int my=(ay+blocking_[1]-1)/blocking_[1];
  int ix,iy,iz;
  block->index().array(&ix,&iy,&iz);
  ix /= blocking_[0];
  iy /= blocking_[1];
  iz /= blocking_[2];
  return (ix + mx*(iy + my*iz));
}

//----------------------------------------------------------------------

DataMsg * MethodOutput::create_data_msg_ (Block * block)
{
  // set face for entire block data
  int if3[3] = {0,0,0};
  // [set child index even though not accessed]
  int ic3[3] = {0,0,0};
  // include ghost zones (all or nothing)
  int g3[3] = {0};
  if (ghost_depth_ > 0) {
    cello::field_descr()->ghost_depth(0,g3,g3+1,g3+2);
  }

  // Create refresh object required by FieldFace
  Refresh * refresh = new Refresh;

  // Initialize refresh fields
  bool any_fields = true;
  if (all_fields_) {
    // all fields
    refresh->add_all_fields();
  } else if (field_list_.size() > 0) {
    // some fields
    refresh->set_field_list(field_list_);
  } else {
    // no fields
    any_fields = false;
  }

  // Initialize refresh particles
  bool any_particles = true;
  if (all_particles_) {
    // all particle types
    refresh->add_all_particles();
  } else if (particle_list_.size() > 0) {
    // some particle types
    refresh->set_particle_list(particle_list_);
  } else {
    // no particle types
    any_particles = false;
  }

  // Create FieldFace object specifying fields to send
  FieldFace * field_face = block->create_face
    (if3,ic3,g3, refresh_same, refresh);

  int gx=-1,gy=-1,gz=-1;
  field_face->ghost(&gx,&gy,&gz);

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
