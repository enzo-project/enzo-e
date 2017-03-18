// See LICENSE_CELLO file for license and copyright information

/// @file     io_Input.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Input class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

Input::Input (const Factory * factory,
	      const FieldDescr * field_descr,
	      const ParticleDescr * particle_descr) throw()
  : file_(0),           // Initialization deferred
    process_(0),        // initialization below
    sync_(0),
    index_charm_(0),
    cycle_(0),
    time_(0),
    file_name_(""),     // set_filename()
    file_args_(),       // set_filename()
    io_block_(0),
    it_field_index_(0),        // set_it_index_field()
    io_field_data_(0),
    it_particle_index_(0),        // set_it_index_particle()
    io_particle_data_(0),
    process_stride_(1) // default one file per process
{

  process_  = CkMyPe();

  io_block_         = factory->create_io_block();
  io_field_data_    = factory->create_io_field_data(field_descr);
  io_particle_data_ = factory->create_io_particle_data(particle_descr);
}

//----------------------------------------------------------------------

Input::~Input () throw()
{
  delete file_;      file_ = 0;
  delete it_field_index_;  it_field_index_ = 0;
  delete it_particle_index_;  it_particle_index_ = 0;
  delete io_block_;  io_block_ = 0;
}

//----------------------------------------------------------------------

void Input::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  PUP::able::pup(p);

  bool up = p.isUnpacking();

  WARNING("Input::pup","skipping file_");
  //    p | *file_;
  p | process_;
  p | sync_;
  p | index_charm_;
  p | cycle_;
  p | time_;
  p | file_name_;
  p | file_args_;
  p | it_field_index_; // PUP::able
  p | it_particle_index_; // PUP::able
  if (up) io_block_ = new IoBlock;
  p | *io_block_;

  WARNING("Input::pup","skipping io_field_data_");
  // if (up) io_field_data_ = new IoFieldData;
  // p | *io_field_data_;

  WARNING("Input::pup","skipping io_particle_data_");
  // if (up) io_particle_data_ = new IoParticleData;
  // p | *io_particle_data_;

  p | process_stride_;
}

//----------------------------------------------------------------------

void Input::set_filename (std::string filename,
			   std::vector<std::string> fileargs) throw()
{
  file_name_ = filename; 
  file_args_ = fileargs;
}

//----------------------------------------------------------------------

void Input::read_simulation
(
 Simulation * simulation
 ) throw()
{
  WARNING("Input::read_simulation()",
	  "This function should not be called");
}

//----------------------------------------------------------------------

void Input::read_hierarchy
(
 Hierarchy * hierarchy,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr
) throw()
{

  if (CkMyPe() == 0)

    // --------------------------------------------------
    // ENTRY: #1 Input::read_hierarchy()-> Block::p_output_read()
    // ENTRY: Block array if Simulation is root
    // --------------------------------------------------
    hierarchy->block_array()->p_output_read (index_charm_);
    // --------------------------------------------------

}

//----------------------------------------------------------------------

Block * Input::read_block
(
 Block * block,
 std::string block_name,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr
 ) throw()
{
  // Read fields

  ItIndex * it_f = it_field_index_;
  if (it_f) {
    for (it_f->first(); ! it_f->done(); it_f->next()  ) {
      int index_field = it_f->value();
      read_field (block->data()->field_data(), 
		  field_descr, index_field);
    }
  }

  ItIndex * it_p = it_particle_index_;
  if (it_p) {
    for (it_p->first(); ! it_p->done(); it_p->next()  ) {
      int index_particle = it_p->value();
      read_particle (block->data()->particle_data(), 
		     particle_descr, index_particle);
    }
  }

  return block;
}

//======================================================================

std::string Input::expand_name_
(
 const std::string              * name_p,
 const std::vector<std::string> * args_p
) const throw()
{
  const std::string & name = *name_p;
  const std::vector<std::string> & args = *args_p;
  
  const int MAX_BUFFER = 255;

  char buffer[MAX_BUFFER+1];
  char buffer_new[MAX_BUFFER+1];

  // Error check no \% in file name

  ASSERT1 ("Input::expand_name_",
	   "File name %s cannot contain '\\%%'",
	   name.c_str(),
	   name.find("\\%") == std::string::npos);

  // Error check variable count equals format conversion specifier count

  std::string rest = name;
  size_t count = 0;
  size_t pos = 0;
  size_t len;
  while ((pos = rest.find("%")) != std::string::npos) {
    count ++;
    len = rest.size();
    rest = rest.substr(pos+1,len-pos-1);
  }

  ASSERT3 ("Input::expand_name_",
	   "The number of format conversion specifiers %d "
	   "associated with file name %s "
	   "must equal the number of variables %d",
	    count, name.c_str(),args.size(),
	   args.size() == count);

  // loop through args[] from the right and replace 
  // format strings with variable values

  std::string left  = name;
  std::string middle = "";
  std::string right = "";

  for (size_t i=0; i<args.size(); i++) {

    // visit variables from right to left
    const std::string & arg = args[args.size() - i - 1];

    size_t pos = left.rfind("%");
    size_t len = left.size();

    middle = left.substr(pos,len-pos);
    left  = left.substr(0,pos);

    strncpy (buffer, middle.c_str(),MAX_BUFFER);
    
    if      (arg == "cycle") { sprintf (buffer_new,buffer, cycle_); }
    else if (arg == "time")  { sprintf (buffer_new,buffer, time_); }
    else if (arg == "proc")  { sprintf (buffer_new,buffer, process_); }
    else 
      {
	ERROR3("Input::expand_name_",
	       "Unknown file variable #%d '%s' for file '%s'",
	       int(i),arg.c_str(),name.c_str());
      }

    right = std::string(buffer_new) + right;

  }

  return left + right;
}

//----------------------------------------------------------------------

void Input::read_meta_ ( meta_type type_meta, Io * io ) throw ()
{
  for (size_t i=0; i<io->meta_count(); i++) {

    void * buffer;
    std::string name;
    int type_scalar;
    int nx,ny,nz;

    // Get object's ith metadata

    io->meta_value(i,& buffer, &name, &type_scalar, &nx,&ny,&nz);

    // Read object's ith metadata
    if ( type_meta == meta_type_group ) {
      file_->group_read_meta(buffer,name.c_str(),&type_scalar,&nx,&ny,&nz);
    } else {
      file_->file_read_meta(buffer,name.c_str(),&type_scalar,&nx,&ny,&nz);
    }
  }
}

