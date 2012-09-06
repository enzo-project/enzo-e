// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputRestart.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-02-02
/// @brief    Implementation of writing restart dumps
///
/// 

#include "io.hpp"

#ifdef CONFIG_USE_CHARM
extern CProxy_SimulationCharm  proxy_simulation;
#endif

//----------------------------------------------------------------------

OutputRestart::OutputRestart
(
 const Factory * factory,
 Parameters * parameters,
 int process_count
) throw ()
  : Output(factory),
    dir_name_(""),
    dir_args_()
{

  set_process_stride(process_count);

#ifndef CONFIG_USE_CHARM
  WARNING("OutputRestart::OutputRestart",
	  "Restart capability only implemented for Charm++");
#endif

  // ASSUMES "Output : <file_group>" is current parameters group
  // set in Problem::initialize_output()

  //--------------------------------------------------
  // parameter: Output : <file_group> : param
  //--------------------------------------------------
  
  if (parameters->type("dir") == parameter_string) {

    dir_name_ = parameters->value_string("dir","");
     
  } else if (parameters->type("dir") == parameter_list) {

    int list_length = parameters->list_length("dir");

    if (list_length > 0) {
      dir_name_ = parameters->list_value_string(0,"dir","");
    }

    for (int index = 1; index<list_length; index++) {
      dir_args_.push_back(parameters->list_value_string(index,"dir",""));
    }

  } else {

      ERROR1("OutputRestart::OutputRestart",
	     "Bad type %d for 'Output : <file_group> : dir' parameter",
	     parameters->type("dir"));

  }

}

//======================================================================

void OutputRestart::write ( const Simulation * simulation ) throw()
{

  TRACE("OutputRestart::write");

#ifdef CONFIG_USE_CHARM

  // Write parameter file

  bool is_root = simulation->group_process()->is_root();

  if (is_root) {

    std::string dir_name = expand_file_name_(&dir_name_,&dir_args_);

    ERROR("OutputRestart::write",
	  "Restart not debugged yet--hangs");
    CkCallback callback(CkIndex_SimulationCharm::s_write(),proxy_simulation);
    CkStartCheckpoint (dir_name.c_str(),callback);
  }
#else

  ERROR("OutputRestart::OutputRestart",
	"Restart capability only implemented for Charm++");
#endif


}

//======================================================================
