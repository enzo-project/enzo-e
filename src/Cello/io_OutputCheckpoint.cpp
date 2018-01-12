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
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 Config * config,
 int process_count
) throw ()
  : Output(index,factory,field_descr,particle_descr),
    restart_file_("")
{

  set_stride_write (process_count);
  
  stride_wait_ = 1;

  TRACE1 ("index = %d",index_);
  TRACE2 ("config->output_dir[%d]=%p",index_,&config->output_dir[index_]);
  TRACE2 ("config->output_dir[%d][0]=%s",
	  index_,config->output_dir[index_][0].c_str());

  restart_file_ = config->restart_file;

}


//----------------------------------------------------------------------

void OutputCheckpoint::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change

  Output::pup(p);

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

  // Testing:time_final

  if (p.type("Testing:time_final") == parameter_list) {
    int length = p.list_length("Testing:time_final");
    config->testing_time_final.resize(length);
    for (int i=0; i<length; i++) {
      config->testing_time_final[i] = p.list_value_float (i,"Testing:time_final",0.0);
    }
  } else {
    config->testing_time_final.resize(1);
    config->testing_time_final[0] = p.value_float
      ("Testing:time_final",0.0);
  }

}

//----------------------------------------------------------------------

void OutputCheckpoint::write_simulation ( const Simulation * simulation ) throw()
{
  TRACE("OutputCheckpoint::write_simulation()");

  std::string dir_name = expand_name_(&dir_name_,&dir_args_);

  simulation->set_phase (phase_restart);

  proxy_main.p_checkpoint(CkNumPes(),dir_name);

}

//======================================================================
