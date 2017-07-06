// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
/// @brief    Implementation of the OutputData class


#include "cello.hpp"
#include "main.hpp"
#include "io.hpp"

//----------------------------------------------------------------------

OutputData::OutputData
(
 int index,
 const Factory * factory,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 Config * config
) throw ()
  : Output(index,factory,field_descr,particle_descr),
    text_block_count_(0)
{
  // Set process stride, with default = 1

  int stride = config->output_stride[index_];

  process_stride_ = stride == 0 ? 1 : stride;
  set_process_stride(process_stride_);

}

//----------------------------------------------------------------------

OutputData::~OutputData() throw()
{
  close();
}

//----------------------------------------------------------------------

void OutputData::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change

  Output::pup(p);

  p | text_block_count_;
}

//======================================================================

void OutputData::open () throw()
{
  std::string file_name = expand_name_(&file_name_,&file_args_);


  close();

  std::string dir = directory();

  Monitor::instance()->print 
    ("Output","writing data file %s",
     (dir + "/" + file_name).c_str());

  file_ = new FileHdf5 (dir,file_name);

  file_->file_create();
}

//----------------------------------------------------------------------

void OutputData::close () throw()
{
  if (file_) file_->file_close();
  delete file_;  file_ = 0;
}

//----------------------------------------------------------------------

void OutputData::finalize () throw ()
{  Output::finalize(); }

//----------------------------------------------------------------------

void OutputData::write_hierarchy
(
 const Hierarchy  * hierarchy,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr
 ) throw()
{
  IoHierarchy io_hierarchy(hierarchy);

  write_meta (&io_hierarchy);

  Output::write_hierarchy(hierarchy, field_descr, particle_descr);

}

//----------------------------------------------------------------------

void OutputData::write_block
( 
  const Block * block,
  const FieldDescr * field_descr,
  const ParticleDescr * particle_descr) throw()
{

  char file[256];
  char dir[256];
  char line[256];

  std::string name_dir = expand_name_(&dir_name_,&dir_args_);

  // Write blocks text file
  
  if (name_dir != "") {

    std::string name_file  = expand_name_(&file_name_,&file_args_);
    const int num_blocks = block->simulation()->hierarchy()->num_blocks();
    int count = 0;
    
    // Write DIR.parameters file
    
    std::string param_file_name = name_dir+"/"+name_dir+".parameters";
    g_parameters.write(param_file_name.c_str(),false);
    
    // Contribute to DIR.block_list file
    
    count = (text_block_count_ == 0) ? num_blocks : 0;
    
    sprintf (file,"%s.block_list",name_dir.c_str());
    sprintf (dir, "%s",           name_dir.c_str());
    sprintf (line,"%s %s\n",      block->name().c_str(),name_file.c_str());
    
    proxy_main.p_text_file_write(strlen(dir)+1, dir,
				 strlen(file)+1, file,
				 strlen(line)+1,line,
				 count);
    
    // Contribute to DIR.file_list file

    if (text_block_count_ == 0) {
      
      count = 0;
    
      sprintf (file,"%s.file_list",name_dir.c_str());
      sprintf (dir, "%s",          name_dir.c_str());
      sprintf (line,"%s\n",        name_file.c_str());
    
      proxy_main.p_text_file_write(strlen(dir)+1, dir,
				   strlen(file)+1, file,
				   strlen(line)+1,line,
				   count);
    }    

    // Increment block counter

    text_block_count_ = (text_block_count_ + 1) % num_blocks;

  }

  // Create file group for block

  std::string group_name = "/" + block->name();

  DEBUG1 ("block name = %s",group_name.c_str());
  file_->group_chdir(group_name);
  file_->group_create();

  // Write block meta data

  io_block()->set_block((Block *)block);

  write_meta_group (io_block());

  // Call write(block) on base Output object

  Output::write_block(block,field_descr,particle_descr);

  file_->group_close();

}

//----------------------------------------------------------------------

void OutputData::write_field_data
( 
  const FieldData * field_data,
  const FieldDescr * field_descr,
  int index_field) throw()
{
  io_field_data()->set_field_descr((FieldDescr*)field_descr);
  io_field_data()->set_field_data((FieldData*)field_data);
  io_field_data()->set_field_index(index_field);

  for (size_t i=0; i<io_field_data()->data_count(); i++) {

    void * buffer;
    std::string name;
    int type;
    int nxd,nyd,nzd;  // Array dimension
    int nx,ny,nz;     // Array size

    // Get ith FieldData data
    io_field_data()->field_array(i, &buffer, &name, &type, 
				 &nxd,&nyd,&nzd,
				 &nx, &ny, &nz);

    // Write ith FieldData data

    file_->mem_create(nx,ny,nz,nx,ny,nz,0,0,0);
    if (nzd > 1) {
      file_->data_create(name.c_str(),type,nzd,nyd,nxd,1,nz,ny,nx,1);
    } else if (nyd > 1) {
      file_->data_create(name.c_str(),type,nyd,nxd,  1,1,ny,nx, 1,1);
    } else {
      file_->data_create(name.c_str(),type,nxd,  1,  1,1,nx,  1,1,1);
    }
    file_->data_write(buffer);
    file_->data_close();
  }

}

//----------------------------------------------------------------------

void OutputData::write_particle_data
( 
  const ParticleData * particle_data,
  const ParticleDescr * particle_descr,
  int it) throw()
{
  const Particle particle ( (ParticleDescr*) particle_descr,
			    (ParticleData*)  particle_data);

  // Write particle data for particle type it

  io_particle_data()->set_particle_descr ( (ParticleDescr*) particle_descr);
  io_particle_data()->set_particle_data  ( (ParticleData*)  particle_data);
  io_particle_data()->set_particle_index(it);

  // loop through attributes 
  // loop through blocks
  // write particle data for [it][ib][ia]

  const int nb = particle.num_batches(it);
  const int na = particle.num_attributes(it);

  // For each particle attribute
  for (int ia=0; ia<na; ia++) {

    int np = particle.num_particles (it);

    const std::string name = "particle "
      +                particle.type_name(it) + " "
      +                particle.attribute_name(it,ia);
    
    const int type = particle.attribute_type(it,ia);

    // create the disk array
    file_->data_create(name.c_str(),type,np,1,1,1,np,1,1,1);
    
    const int dx = particle.stride(it,ia);
    int i0 = 0;

    // for each batch of particles
    
    for (int ib=0; ib<nb; ib++) {
      
      int mb = particle.num_particles(it,ib);

      // create the memory space for the batch
      file_->mem_create(mb,1,1,mb,1,1,0,0,0);
      
      const void * buffer = (const void *) particle.attribute_array(it,ia,ib);

      // find the hyper_slab of the disk dataset
      file_->data_slice
	(np, 1, 1, 1,
	 mb, 1, 1, 1,
	 i0, 0, 0, 0);
      
      i0 += mb;

      // write the batch to disk
      file_->data_write(buffer);
      
      if (mb > 0) file_->mem_close();
    }

    // check that the number of particles equals the number written
    
    ASSERT2 ("OutputData::write_particle_data()",
	     "Particle count mismatch %d particles %d written",
	     np,i0,
	     np == i0);

    // close the attribute dataset
    if (np > 0) file_->data_close();
  }

}

//======================================================================
