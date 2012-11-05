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
 int index,
 const Factory * factory,
 Config * config,
 int process_count
) throw ()
  : Output(index,factory),
    dir_name_(""),
    dir_args_()
{

  set_process_stride(process_count);

#ifndef CONFIG_USE_CHARM
  WARNING("OutputRestart::OutputRestart",
	  "Restart capability only implemented for Charm++");
#endif

  TRACE1 ("index = %d",index_);
  TRACE2 ("config->output_dir[%d]=%p",index_,&config->output_dir[index_]);
  TRACE2 ("config->output_dir[%d][0]=%s",
	  index_,config->output_dir[index_][0].c_str());
  dir_name_ = config->output_dir[index_][0];
  TRACE0;

  for (size_t i=1; i<config->output_dir[index_].size(); i++) {
    dir_args_.push_back(config->output_dir[index_][i]);
  }

}


//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void OutputRestart::pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    Output::pup(p);

    p | dir_name_;
    p | dir_args_;
  }

#endif

//======================================================================

void OutputRestart::write_simulation ( const Simulation * simulation ) throw()
{

  TRACE("OutputRestart::write_simulation()");

#ifdef CONFIG_USE_CHARM

  // Write parameter file

  bool is_root = simulation->group_process()->is_root();

  if (is_root) {

    std::string dir_name = expand_file_name_(&dir_name_,&dir_args_);
    char dir_char[255];
    strcpy(dir_char,dir_name.c_str());
    // ERROR("OutputRestart::write_simulation",
    // 	  "Restart not debugged yet--hangs");
    CkCallback callback(CkIndex_SimulationCharm::s_write(),proxy_simulation);
    
    CkStartCheckpoint (dir_char,callback);
  }
#else

  ERROR("OutputRestart::OutputRestart",
	"Restart capability only implemented for Charm++");
#endif


}

//======================================================================
