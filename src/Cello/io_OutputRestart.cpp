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
  : OutputData(factory),
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
  bool is_restart = type == "restart";
  int init_cycle  = parameters->value_integer ("Initial:cycle",-1);

  TRACE3("type = %s is_restart %d init_cycle %d",
	 type.c_str(),is_restart,init_cycle);
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

void OutputRestart::write_simulation
( const Simulation * simulation ) throw()
/// Note factory, field_descr, and hierarchy needed since otherwise
/// Simulation functions must be called, which would introduce a circular
/// dependence between Simulation and Output components
{

  // Write parameter file

  bool is_root = simulation->group_process()->is_root();

  if (is_root) {

    std::string param_name = expand_file_name(&param_name_,&param_args_);

    Parameters * parameters = simulation->parameters();

    // Update Initial parameters

  //--------------------------------------------------
  // parameter: Input : type
  // parameter: Input : name
  // parameter: Input : param
  // parameter: Input : cycle
  // parameter: Input : time
  //--------------------------------------------------

    parameters->set_string  ("Initial:type", "restart");
    WARNING("OutputRestart::write_simulation",
	    "Initial:name setting requires set_list() capability");
    parameters->set_string  ("Initial:name", file_name_.c_str());

    parameters->set_integer ("Initial:cycle",simulation->cycle());
    parameters->set_float   ("Initial:time", simulation->time());

    // Write restart parameter file
    parameters->write(param_name.c_str());
    
  }

  Output::write_simulation(simulation);

}

//======================================================================
