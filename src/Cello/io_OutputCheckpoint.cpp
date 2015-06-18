// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputCheckpoint.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-02-02
/// @brief    Implementation of writing checkpoint dumps
   
#include "io.hpp"
#include "main.hpp"

extern CProxy_Simulation  proxy_simulation;

//----------------------------------------------------------------------

OutputCheckpoint::OutputCheckpoint
(
 int index,
 const Factory * factory,
 Config * config,
 int process_count
) throw ()
  : Output(index,factory),
    dir_name_(""),
    dir_args_(),
    restart_file_("")
{

  set_process_stride(process_count);

  TRACE1 ("index = %d",index_);
  TRACE2 ("config->output_dir[%d]=%p",index_,&config->output_dir[index_]);
  TRACE2 ("config->output_dir[%d][0]=%s",
	  index_,config->output_dir[index_][0].c_str());
  dir_name_ = config->output_dir[index_][0];
  TRACE0;

  for (size_t i=1; i<config->output_dir[index_].size(); i++) {
    dir_args_.push_back(config->output_dir[index_][i]);
  }

  restart_file_ = config->restart_file;

}


//----------------------------------------------------------------------

void OutputCheckpoint::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change

  Output::pup(p);

  p | dir_name_;
  p | dir_args_;
  p | restart_file_;

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  const bool l_unpacking = p.isUnpacking();
  const bool l_restarting = simulation && simulation->phase()==phase_restart;
  const bool l_restart_file = (restart_file_ != "");
  if (l_unpacking && l_restarting && l_restart_file) {
    update_config_();
  }
}

//======================================================================

void OutputCheckpoint::update_config_()
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  Parameters p (simulation->monitor());
  p.read(restart_file_.c_str());

  Config * config = (Config *) simulation->config();

  // Field:courant

  config->field_courant = p.value_float
    ("Field:courant",config->field_courant);

  // Testing:time_final

  config->testing_time_final = p.value_float
    ("Testing:time_final",config->testing_time_final);

}

//----------------------------------------------------------------------

void OutputCheckpoint::write_simulation ( const Simulation * simulation ) throw()
{
  TRACE("OutputCheckpoint::write_simulation()");

  std::string dir_name = expand_file_name_(&dir_name_,&dir_args_);

  simulation->set_phase (phase_restart);

  proxy_main.p_checkpoint(CkNumPes(),dir_name);

}

//======================================================================
