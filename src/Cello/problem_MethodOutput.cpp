// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOutput.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-04-27
/// @brief    Implementation of the Output "method"

#include "problem.hpp"
#include "charm_simulation.hpp"
#include "test.hpp"

#define TRACE_METHOD_OUTPUT
#define TRACE_FILE

//--------------------------------------------------

#ifdef TRACE_METHOD_OUTPUT
#   undef TRACE_METHOD_OUTPUT
#   define TRACE_METHOD_OUTPUT(block,msg,tag)                              \
  {                                                                     \
    CkPrintf ("TRACE_METHOD_OUTPUT %s %s %s\n",std::string(msg).c_str(),block->name().c_str(),tag); \
    fflush(stdout);                                                     \
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

//----------------------------------------------------------------------

MethodOutput::MethodOutput
(std::vector< std::string > file_name,
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
    is_count_(-1)
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

void MethodOutput::compute ( Block * block) throw()
{
  TRACE_METHOD_OUTPUT(block,"compute","...");

  // barrier to ensure tree traversal doesn't reach a block
  // before the method has started

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
    msg_output->set_index_send (block->index());
    // Get its block_trace object
    BlockTrace * block_trace = msg_output->block_trace();

    TRACE_FILE(block,file,"write meta hierarchy");
    file_write_hierarchy_(file);
    
    if (block->is_leaf()) {

      TRACE_FILE(block,file,std::string("WRITE ")+block->name().c_str());
    }

    
    if ( ! block_trace->next(block->is_leaf())) {

      Index index_next = block_trace->top();
      TRACE_METHOD_OUTPUT(block,"continue next",msg_output->tag());
      cello::block_array()[index_next].p_method_output_next(msg_output);
      
    } else {
      // Close file
      FileHdf5 * file = msg_output->file();
      file->file_close();
      delete file;
      // Delete message
      
      TRACE_FILE(block,file,"close");
      TRACE_METHOD_OUTPUT(block,"done",msg_output->tag());

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
  msg_output->method_output_->next(this,msg_output);
}

//----------------------------------------------------------------------

void MethodOutput::next(Block * block, MsgOutput * msg_output_in )
{
  // Copy incoming message and delete old (cannot reuse!)
  MsgOutput * msg_output = new MsgOutput (*msg_output_in);
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
    msg_output->set_block(block);
    cello::block_array()[index_home].p_method_output_write(msg_output);
  } else {
    // if non-leaf, forward to child
    Index index_next = bt->top();
    TRACE_METHOD_OUTPUT(block,"next next",msg_output->tag());
    cello::block_array()[index_next].p_method_output_next(msg_output);
  }
  block->compute_done();

}

//----------------------------------------------------------------------

void Block::p_method_output_write (MsgOutput * msg_output)
{
  msg_output->method_output_->write(this,msg_output);
}

//----------------------------------------------------------------------

void MethodOutput::write(Block * block, MsgOutput * msg_output_in )
{
  // Copy incoming message and delete old (cannot reuse!)
  MsgOutput * msg_output = new MsgOutput (*msg_output_in);
  delete msg_output_in;
  msg_output_in = nullptr;

  TRACE_METHOD_OUTPUT(block,"write",msg_output->tag());
  // delete data
  FileHdf5 * file = msg_output->file();

  file_write_block_(file,block,msg_output);

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
  int ax,ay,az;
  cello::hierarchy()->root_blocks(&ax,&ay,&az);
  int mx=(ax+blocking_[0]-1)/blocking_[0];
  int my=(ay+blocking_[1]-1)/blocking_[1];
  int ix,iy,iz;
  block->index().array(&ix,&iy,&iz);
  ix /= blocking_[0];
  iy /= blocking_[1];
  iz /= blocking_[2];
  int count = ix + mx*(iy + my*iz);
  std::string file_name = cello::expand_name(&file_name_,count,block);

  // Create File
  FileHdf5 * file = new FileHdf5 (path_name, file_name);
  TRACE_FILE(block,file,std::string("open ")+path_name+" "+file_name);
  file->file_create();
  return file;
}

//----------------------------------------------------------------------

void MethodOutput::file_write_hierarchy_(FileHdf5 * file)
{
  IoHierarchy io_hierarchy(cello::hierarchy());
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
  TRACE_FILE(block,file,"write ");
  ASSERT2("MethodOutput::file_write_block",
          "Block mismatch %s != %s",
          block->name(msg_output->index_send()).c_str(),
          msg_output->block_name.c_str(),
          (block->name(msg_output->index_send()) == 
           msg_output->block_name));
          
}
