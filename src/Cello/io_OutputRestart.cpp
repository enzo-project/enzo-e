// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputRestart.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-02-02
/// @brief    Implementation of writing restart dumps
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

OutputRestart::OutputRestart
(
 const Factory * factory,
 Parameters * parameters
) throw ()
  : OutputData(factory,parameters),
    param_name_(""),
    param_args_()
{

  // ASSUMES "Output : <file_group>" is current parameters group
  // set in Problem::initialize_output()

  //--------------------------------------------------
  // parameter: Output : <file_group> : param
  //--------------------------------------------------

  
  if (parameters->type("param") == parameter_string) {

    param_name_ = parameters->value_string("param","");
     
  } else if (parameters->type("param") == parameter_list) {

    int list_length = parameters->list_length("param");

    if (list_length > 0) {
      param_name_ = parameters->list_value_string(0,"param","");
    }

    for (int index = 1; index<list_length; index++) {
      param_args_.push_back(parameters->list_value_string(index,"param",""));
    }

  } else {

      ERROR1("OutputRestart::OutputRestart",
	     "Bad type %d for 'Output : <file_group> : param' parameter",
	     parameters->type("param"));

  }

  //--------------------------------------------------
  // parameter: Initial : type
  // parameter: Initial : cycle
  //--------------------------------------------------

  // Skip first cycle for restart if this is a restart

  std::string type = parameters->value_string("Initial:type","");
  bool is_restart = (type == "restart");
  int init_cycle  = parameters->value_integer ("Initial:cycle",-1);

  if (is_restart) {
    schedule()->set_skip_cycle(init_cycle);
  }

}

//======================================================================

void OutputRestart::finalize () throw ()
{
  OutputData::finalize();
}

//----------------------------------------------------------------------

void OutputRestart::write ( const Simulation * simulation ) throw()
/// Note factory, field_descr, and hierarchy needed since otherwise
/// Simulation functions must be called, which would introduce a circular
/// dependence between Simulation and Output components
{

  // Write parameter file

  bool is_root = simulation->group_process()->is_root();

  if (is_root) {

    std::string param_name = expand_file_name_(&param_name_,&param_args_);

    Parameters * parameters = simulation->parameters();

    // Update Initial parameters

  //--------------------------------------------------
  // parameter: Initial : type
  // parameter: Initial : name
  // parameter: Initial : param
  // parameter: Initial : cycle
  // parameter: Initial : time
  //--------------------------------------------------

    parameters->set_string  ("Initial:type", "restart");
    
    int list_length = 1 + file_args_.size();
    parameters->set_list_length  ("Initial:name", list_length);
    parameters->set_list_string (0,"Initial:name",file_name_.c_str());
    for (size_t i=0; i<file_args_.size(); i++) {
      parameters->set_list_string (i+1,"Initial:name",file_args_[i].c_str());
    }
    parameters->done_set_list("Initial:name");

    parameters->set_integer ("Initial:cycle",simulation->cycle());
    parameters->set_float   ("Initial:time", simulation->time());

    // Write restart parameter file
    parameters->write(param_name.c_str());
    
  }

  // Call base write()
  write_ (simulation);

}

//======================================================================
