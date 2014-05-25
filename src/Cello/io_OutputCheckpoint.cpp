// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputCheckpoint.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-02-02
/// @brief    Implementation of writing checkpoint dumps
///
/// 


   
#include "io.hpp"
#include "main.hpp"

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
    dir_args_()
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

}


//----------------------------------------------------------------------

void OutputCheckpoint::pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    Output::pup(p);

    p | dir_name_;
    p | dir_args_;
  }

//======================================================================

void OutputCheckpoint::write_simulation ( const Simulation * simulation ) throw()
{
  TRACE("OutputCheckpoint::write_simulation()");

  std::string dir_name = expand_file_name_(&dir_name_,&dir_args_);

  proxy_main.p_checkpoint(CkNumPes(),dir_name);

}

//======================================================================
