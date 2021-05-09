// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOutput.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-04-27
/// @brief    Implementation of the Output "method"

#include "problem.hpp"
#include "charm_simulation.hpp"
#include "test.hpp"

// #define TRACE_METHOD_OUTPUT
// #define TRACE_FILE
// #define TRACE_MSG_OUTPUT
// #define TRACE_FIELD_SUMS /* use for cello-bug/210503 checking */

#define TRACE_BLOCK_DENSITY(BLOCK_NAME,DATA)                           \
  {                                                                     \
    Field field = DATA->field();                                        \
    if (BLOCK_NAME == "B000_000") {                                     \
      int index_density = field.field_id("density");                    \
      int mx,my,mz;                                                     \
      field.dimensions(index_density,&mx,&my,&mz);                      \
      cello_float * density = (cello_float *)field.values(index_density); \
      CkPrintf ("TRACE_BLOCK_DENSITY %s %d block %s density [4,4] %p %g\n", \
        __FILE__,__LINE__,BLOCK_NAME.c_str(), (void *)&(density[4+mx*(4)]), density[4+mx*(4)]); \
    }                                                                   \
  }

//--------------------------------------------------

#ifdef TRACE_METHOD_OUTPUT
#   undef TRACE_METHOD_OUTPUT
#   define TRACE_METHOD_OUTPUT(block,msg,tag)                           \
  {                                                                     \
  CkPrintf ("TRACE_METHOD_OUTPUT %s %s %s\n",std::string(msg).c_str(),block->name().c_str(),tag); \
  fflush(stdout);                                                       \
}
#else
#   define TRACE_METHOD_OUTPUT(block,msg,tag) /* ... */
#endif

//--------------------------------------------------

#ifdef TRACE_FILE
#   undef TRACE_FILE
#   define TRACE_FILE(block,file,msg)           \
  {                                             \
    CkPrintf ("TRACE_FILE %p %s %s\n",          \
              (void *)file,                     \
              block->name().c_str(),            \
              std::string(msg).c_str());        \
    fflush(stdout);                             \
  }
#else
#   define TRACE_FILE(block,file,msg) /* ... */
#endif

//--------------------------------------------------

#ifdef TRACE_MSG_OUTPUT
#   undef TRACE_MSG_OUTPUT                                            
#   define TRACE_MSG_OUTPUT(MSG_OUTPUT,MSG)                     \
  CkPrintf ("TRACE_MSG_OUTPUT %s:%d %d %p %s %s\n",                \
            __FILE__,__LINE__,CkMyPe(),(void *) MSG_OUTPUT,MSG_OUTPUT->tag(),MSG);  
#else  
#   define TRACE_MSG_OUTPUT(MSG_OUTPUT,MSG) /* ... */
#endif
    

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
   int blocking_x,
   int blocking_y,
   int blocking_z)
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
      factory_(factory)
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

  if (ghost_depth > 0) refresh->set_ghost_depth (ghost_depth);
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

  ASSERT("MethodOutput()",
         "MethodOutput requires either field_list or particle_list "
         "to be non-empty",
         (all_fields || all_particles ||
          field_list.size() > 0 || particle_list.size() > 0));
         
  // initialize blocking partition
  blocking_[0] = blocking_x;
  blocking_[1] = blocking_y;
  blocking_[2] = blocking_z;

  is_count_ = cello::scalar_descr_int()->new_value("method_output:count");

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
  p | factory_nonconst_;

}

//----------------------------------------------------------------------

void MethodOutput::compute ( Block * block) throw()
{
  TRACE_METHOD_OUTPUT(block,"compute","...");

  // barrier to ensure tree traversal doesn't reach a block
  // before the method has started

  TRACE_BLOCK_DENSITY(block->name(),block->data());
  
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
  TRACE_METHOD_OUTPUT(block,"compute_continue","...");

  if (is_writer_(block->index())) {
    int a3[3];
    block->index().array(a3,a3+1,a3+2);
    int index_min[3] = {a3[0] - (a3[0] % blocking_[0]),
                        a3[1] - (a3[1] % blocking_[1]),
                        a3[2] - (a3[2] % blocking_[2])};
    int index_max[3] = {index_min[0] + blocking_[0],
                        index_min[1] + blocking_[1],
                        index_min[2] + blocking_[2]};
    BlockTrace bt (cello::rank(),index_min,index_max);

    // Create and open file
    FileHdf5 * file = file_open_(block,a3);
    
    // Create output message
    MsgOutput * msg_output = new MsgOutput(bt,this,file);
    TRACE_MSG_OUTPUT(msg_output,"new");
    msg_output->set_index_send (block->index());
    // Get its block_trace object
    BlockTrace * block_trace = msg_output->block_trace();

    file_write_hierarchy_(file);
    
    if (block->is_leaf()) {

      TRACE_FILE(block,file,std::string("WRITE ")+block->name().c_str());
      msg_output->set_block(block,factory_);
      file_write_block_(file,block,nullptr);
      
    }

    
    if ( ! block_trace->next(block->is_leaf())) {
      Index index_next = block_trace->top();
      TRACE_METHOD_OUTPUT(block,"continue next",msg_output->tag());
      msg_output->del_block();
      cello::block_array()[index_next].p_method_output_next(msg_output);
      
    } else {
      // Close file
      FileHdf5 * file = msg_output->file();
      TRACE_FILE(block,file,"close file");
      file->file_close();
      delete file;
      TRACE_METHOD_OUTPUT(block,"done",msg_output->tag());

      // Delete message
      TRACE_MSG_OUTPUT(msg_output,"delete");
      delete msg_output;
      msg_output = nullptr;

      // Exit MethodOutput
      block->compute_done();
    }
  }

  if (block->level() < 0) {
    // negative level blocks do not participate
    TRACE_METHOD_OUTPUT(block,"done","...");
    block->compute_done();
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
{
  // Copy incoming message and delete old (cannot reuse!)
  MsgOutput * msg_output = new MsgOutput (*msg_output_in);
  TRACE_MSG_OUTPUT(msg_output,"new");

  TRACE_MSG_OUTPUT(msg_output_in,"delete");
  delete msg_output_in;
  msg_output_in = nullptr;
  
  TRACE_METHOD_OUTPUT(block,"next",msg_output->tag());
  const bool is_leaf = block->is_leaf();
  BlockTrace * bt = msg_output->block_trace();
  bt->next(block->is_leaf());
  msg_output->set_index_send (block->index());

  if (is_leaf) {
    // if leaf, copy data to msg and send to writer
    Index index_home = bt->home();
    DataMsg * data_msg = create_data_msg_(block);
    msg_output->set_data_msg(data_msg);
    msg_output->set_block(block,factory_);
    msg_output->io_block()->print("next");
    cello::block_array()[index_home].p_method_output_write(msg_output);
  } else {
    // if non-leaf, forward to child
    Index index_next = bt->top();
    TRACE_METHOD_OUTPUT(block,"next next",msg_output->tag());
    msg_output->del_block();
    cello::block_array()[index_next].p_method_output_next(msg_output);
  }
  block->compute_done();

}

//----------------------------------------------------------------------

void MethodOutput::write(Block * block, MsgOutput * msg_output_in )
{
  // Write the block data
  
  FileHdf5 * file = msg_output_in->file();
  TRACE_FILE(block,file,std::string("WRITE ")+block->name().c_str());
  file_write_block_(file,block,msg_output_in);

  // Copy incoming message
  MsgOutput * msg_output = new MsgOutput (*msg_output_in);
  TRACE_MSG_OUTPUT(msg_output,"new");
  // and delete it (cannot reuse)
  TRACE_MSG_OUTPUT(msg_output_in,"delete");
  delete msg_output_in;
  msg_output_in = nullptr;

  msg_output->set_index_send (block->index());
  BlockTrace * bt = msg_output->block_trace();
  Index index_next = bt->top();
  Index index_home = bt->home();
  if (index_next == index_home) {
    // done
    TRACE_FILE(block,file,"close");
    block->compute_done();
  } else {
    TRACE_METHOD_OUTPUT(block,"write next",msg_output->tag());
    msg_output->del_block();
    cello::block_array()[index_next].p_method_output_next(msg_output);
  }  
}

//----------------------------------------------------------------------

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
  std::string path_name = cello::directory(&path_name_,(*counter),block);
  (*counter)++;

  // Generate file name
  int count = file_count_(block);
  std::string file_name = cello::expand_name(&file_name_,count,block);

  // Create File
  FileHdf5 * file = new FileHdf5 (path_name, file_name);
  TRACE_FILE(block,file,std::string("file create ")+path_name+" "+file_name);
  file->file_create();
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
    block->name() : msg_output->block_name;

  IoBlock * io_block = (is_local) ?
    factory_->create_io_block() : msg_output->io_block();
  if (is_local) io_block->set_block(block);

  double block_lower[3];
  double block_upper[3];
  if (is_local) {
    block->lower(block_lower,block_lower+1,block_lower+2);
    block->upper(block_upper,block_upper+1,block_upper+2);
  }
    
  double * lower = (is_local)?
    block_lower : msg_output->block_lower;
  double * upper = (is_local)?
    block_upper : msg_output->block_upper;

  // Generate path name

  const int count = file_count_(block);
  std::string path_name = cello::expand_name(&path_name_,count,block);
  std::string file_name = cello::expand_name(&file_name_,count,block);

  // Create file group for block

  std::string group_name = "/" + block_name;

  file->group_chdir(group_name);
  file->group_create();

  // Write block meta data

  io_block->print("write");
  write_meta_ (file, io_block, "group");

  // Call write(block) on base Output object

  // Create new data object to hold MsgOutput/DataMsg fields and particles
  
  int nx,ny,nz;
  int num_field_data = 1;
  block->data()->field().size(&nx,&ny,&nz);

  Data * data;

  if (is_local) {
    data = block->data();
  } else {
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
  
  FieldData * field_data = data->field_data();

  // Update data fields with fields from MsgOutput

  TRACE_BLOCK_DENSITY (block_name,data);
  // Write the Block data
  for (int i_f=0; i_f<field_list_.size(); i_f++) {
    const int index_field = field_list_[i_f];
    IoFieldData * io_field_data = factory_->create_io_field_data();

    void * buffer;
    std::string name;
    int type;
    int mx,my,mz;  // Array dimension
    int nx,ny,nz;     // Array size

    io_field_data->set_field_data((FieldData*)field_data);
    io_field_data->set_field_index(index_field);
    io_field_data->field_array
      (&buffer, &name, &type, &mx,&my,&mz, &nx,&ny,&nz);

#ifdef TRACE_FIELD_SUMS    
  {
    Field field = block->data()->field();
    int index_field = field.field_id("density");
    int mx,my,mz;
    int gx,gy,gz;
    field.dimensions(index_field,&mx,&my,&mz);
    field.ghost_depth(index_field,&gx,&gy,&gz);
    long double sum = 0;
    cello_float * temp = (cello_float *)
      field_data->values(cello::field_descr(),index_field);
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
        for (int ix=gx; ix<mx-gx; ix++) {
          int i0 = ix + mx*(iy + my*iz);
          CkPrintf ("DEBUG_VALUE RECV %s %d %d %d %g\n",block_name.c_str(),
                    ix,iy,iz,temp[i0]);
          fflush(stdout);
          sum += temp[i0];
        }
      }
    }
    CkPrintf ("DEBUG_WRITE_BLOCK VALUE RECV %s %d field %d sum %Lg comm %d\n",
              block_name.c_str(), CkMyPe(),
              index_field, sum,
              (block->name() == block_name) ? 0 : 1);
  }
#endif 
 
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
  // Output::write_block(block);
  // ItIndex * it_f = it_field_index_;
  // if (it_f) {
  //   for (it_f->first(); ! it_f->done();  it_f->next()  ) {
  //     const FieldData * field_data = msg_output
  //     write_field_data (field_data, it_f->value());
  //   }
  // }

  // // Write particles

  // ItIndex * it_p = it_particle_index_;
  // if (it_p) {
  //   for (it_p->first(); ! it_p->done();  it_p->next()  ) {
  //     const ParticleData * particle_data = block->data()->particle_data();
  //     write_particle_data (particle_data, it_p->value());
  //   }
  // }

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
  // Create FieldFace for interpolating field data to child ic3[]

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

#ifdef TRACE_FIELD_SUMS    
  {
    Field field = block->data()->field();
    int index_field = field.field_id("density");
    int mx,my,mz;
    int gx,gy,gz;
    field.dimensions(index_field,&mx,&my,&mz);
    field.ghost_depth(index_field,&gx,&gy,&gz);
    long double sum = 0;
    cello_float * temp = (cello_float *) field.values(index_field);
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
        for (int ix=gx; ix<mx-gx; ix++) {
          int i0 = ix + mx*(iy + my*iz);
          sum += temp[i0];
          CkPrintf ("DEBUG_VALUE SEND %s %d %d %d %g\n",block->name().c_str(),
                    ix,iy,iz,temp[i0]);
          fflush(stdout);
        }
      }
    }
    CkPrintf ("DEBUG_WRITE_BLOCK VALUE SEND %s %d field %d sum %Lg\n",
              block->name().c_str(), CkMyPe(),
              index_field, sum);
  }
#endif  

  // Create FieldFace object specifying fields to send
  FieldFace * field_face = block->create_face 
    (if3,ic3,g3, refresh_same, refresh, true);

  // Create data message object to send
  DataMsg * data_msg = new DataMsg;
  data_msg -> set_field_face (field_face,true);
  data_msg -> set_field_data (block->data()->field_data(),false);
  return data_msg;
}

