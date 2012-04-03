// See LICENSE_CELLO file for license and copyright information

/// @file     io_Output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Output class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

Output::Output (const Factory * factory) throw()
  : file_(0),           // Initialization deferred
    schedule_(new Schedule),
    process_(0),        // initialization below
#ifdef CONFIG_USE_CHARM
    counter_(1),        // default process-per-stride
    index_charm_(0),
#endif
    cycle_(0),
    count_output_(0),
    time_(0),
    file_name_(""),     // set_filename()
    file_args_(),       // set_filename()
    it_field_(0),        // set_it_field()
    io_block_(0),
    io_field_block_(0),
    process_stride_(1) // default one file per process

{

  GroupProcess * group_process = GroupProcess::create();

  process_  = group_process->rank();
  delete group_process;

  io_block_       = factory->create_io_block();
  io_field_block_ = factory->create_io_field_block();
}

//----------------------------------------------------------------------

Output::~Output () throw()
{
  delete schedule_;  schedule_ = 0;
  delete file_;      file_ = 0;
  delete it_field_;  it_field_ = 0;
  delete io_block_;  io_block_ = 0;
}
//----------------------------------------------------------------------

bool Output::is_scheduled (int cycle, double time)
{
  cycle_ = cycle;
  time_  = time;
  return schedule()->write_this_cycle(cycle_,time_);
}

//----------------------------------------------------------------------

double Output::update_timestep (double time, double dt) const throw ()
{
  return schedule_->update_timestep(time,dt); 
}

//======================================================================

std::string Output::expand_file_name_
(
 const std::string              * file_name_p,
 const std::vector<std::string> * file_args_p
) const throw()
{
  const std::string & file_name = *file_name_p;
  const std::vector<std::string> & file_args = *file_args_p;
  
  const int MAX_BUFFER = 255;

  char buffer[MAX_BUFFER];
  char buffer_new[MAX_BUFFER];

  // Error check no \% in file name

  ASSERT1 ("Output::expand_file_name_",
	   "File name %s cannot contain '\\%%'",
	   file_name.c_str(),
	   file_name.find("\\%") == std::string::npos);

  // Error check variable count equals format conversion specifier count

  std::string file_rest = file_name;
  size_t count = 0;
  size_t pos = 0;
  size_t len;
  while ((pos = file_rest.find("%")) != std::string::npos) {
    count ++;
    len = file_rest.size();
    file_rest = file_rest.substr(pos+1,len-pos-1);
  }

  ASSERT3 ("Output::expand_file_name_",
	   "The number of format conversion specifiers %d "
	   "associated with file name %s "
	   "must equal the number of variables %d",
	    count, file_name.c_str(),file_args.size(),
	   file_args.size() == count);

  // loop through file_args[] from the right and replace 
  // format strings with variable values

  std::string file_left  = file_name;
  std::string file_middle = "";
  std::string file_right = "";

  for (size_t i=0; i<file_args.size(); i++) {

    // visit variables from right to left
    const std::string & arg = file_args[file_args.size() - i - 1];

    size_t pos = file_left.rfind("%");
    size_t len = file_left.size();

    file_middle = file_left.substr(pos,len-pos);
    file_left  = file_left.substr(0,pos);

    strcpy (buffer, file_middle.c_str());
    
    if      (arg == "cycle") { sprintf (buffer_new,buffer, cycle_); }
    else if (arg == "time")  { sprintf (buffer_new,buffer, time_); }
    else if (arg == "count") { sprintf (buffer_new,buffer, count_output_); }
    else if (arg == "proc")  { sprintf (buffer_new,buffer, process_); }
    else 
      {
	ERROR3("Output::expand_file_name_",
	       "Unknown file variable #%d '%s' for file '%s'",
	       int(i),arg.c_str(),file_name.c_str());
      }

    file_right = std::string(buffer_new) + file_right;

  }

  return file_left + file_right;
}

//----------------------------------------------------------------------

void Output::write_meta_ ( meta_type type_meta, Io * io ) throw ()
{
  for (size_t i=0; i<io->meta_count(); i++) {

    void * buffer;
    std::string name;
    scalar_type type_scalar;
    int nx,ny,nz;

    // Get object's ith metadata

    io->meta_value(i,& buffer, &name, &type_scalar, &nx,&ny,&nz);

    // Write object's ith metadata
    if ( type_meta == meta_type_group ) {
      file_->group_write_meta(buffer,name.c_str(),type_scalar,nx,ny,nz);
    } else {
      file_->file_write_meta(buffer,name.c_str(),type_scalar,nx,ny,nz);
    }
  }
}

//----------------------------------------------------------------------

void Output::write_simulation
( const Simulation * simulation ) throw()
{
  write_hierarchy (simulation->hierarchy(), simulation->field_descr());
}

//----------------------------------------------------------------------

void Output::write_hierarchy
( 
 const Hierarchy * hierarchy,
 const FieldDescr * field_descr
  ) throw()
{

  ItPatch it_patch (hierarchy);

  // (*) write data patch_list_

  while (const Patch * patch = ++it_patch) {

    // NO OFFSET: ASSUMES ROOT PATCH
    write_patch (patch, field_descr, 0,0,0);

  }
}

//----------------------------------------------------------------------

void Output::write_patch 
(
 const Patch * patch,
 const FieldDescr * field_descr,
 int ixp0, int iyp0, int izp0
 ) throw()

{

#ifdef CONFIG_USE_CHARM

  // CHARM++ Block callback for write_block()

  if (patch->blocks_allocated()) {
    patch->block_array().entry_write (index_charm_);
  }

#else

  ItBlock it_block (patch);
  while (const Block * block = ++it_block) {

    // NO OFFSET: ASSUMES ROOT PATCH
    write_block (block, field_descr, 0,0,0);

  }
#endif
}

//----------------------------------------------------------------------

void Output::write_block
(
 const Block * block,
 const FieldDescr * field_descr,
 int ixp0, int iyp0, int izp0
 ) throw()
{
  // Write fields

  for (it_field_->first(); ! it_field_->done(); it_field_->next()  ) {
    const FieldBlock * field_block = block->field_block();
    write_field (field_block,  field_descr, it_field_->value());
  }
}

