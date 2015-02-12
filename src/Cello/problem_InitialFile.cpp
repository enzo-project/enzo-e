// See LICENSE_CELLO file for license and copyright information

/// @file     problem_InitialFile.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-02-16
/// @brief    Implementation of the InitialFile class

#include "cello.hpp"

#include "problem.hpp"

//----------------------------------------------------------------------

InitialFile::InitialFile
(Parameters * parameters,
 int cycle, double time) throw ()
  : Initial (cycle,time),
    parameters_(parameters),
    input_(0),
    block_sync_(0)
{
}

//----------------------------------------------------------------------

InitialFile::~InitialFile() throw()
{
  delete input_; input_ = 0;
}

//----------------------------------------------------------------------

void InitialFile::pup (PUP::er &p)
{
  TRACEPUP;

  Initial::pup(p);

  bool up = p.isUnpacking();

  if (up) parameters_ = new Parameters;
  p | *parameters_;

  p | input_; // PUP::able

  p | block_sync_;

}

//----------------------------------------------------------------------

void InitialFile::enforce_block
(
 Block            * block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()
{
  ASSERT ("InitialFile::enforce_block",
	  "Input block is expected to be NULL",
	  block == 0);

  if (! input_) input_ = new InputData(hierarchy->factory());

  if (! input_->is_open() ) {

    std::string              file_name = "";
    std::vector<std::string> file_args;

    get_filename_(&file_name,&file_args);

    input_->open();

  }

  // input_->read_hierarchy(hierarchy,field_descr);

  INCOMPLETE("InitialFile::enforce_block");

  File * file = input_->file();
  int num_blocks = file->group_count();
  block_sync_.set_stop(num_blocks);

  for (int i = 0; i<num_blocks; i++) {

    std::string block_name = file->group_name(i).c_str();

    input_->read_block(0,block_name,field_descr);

  }
}

//----------------------------------------------------------------------

void InitialFile::get_filename_
(
 std::string * file_name,
 std::vector<std::string> * file_args
 ) throw()
{
  // parameter: Initial : name

  parameters_->group_set(0,"Initial");

  if (parameters_->type("name") == parameter_string) {

    *file_name = parameters_->value_string("name","");
     
  } else if (parameters_->type("name") == parameter_list) {

    int list_length = parameters_->list_length("name");

    *file_name = parameters_->list_value_string(0,"name","");

    for (int index = 1; index<list_length; index++) {
      file_args->push_back(parameters_->list_value_string(index,"name",""));
    }

  } else {

    ERROR1("InitialFile::enforce_block",
	   "Bad type %d for 'Initial : name' parameter",
	   parameters_->type("name"));

  }

  input_->set_filename (*file_name, *file_args);

}
