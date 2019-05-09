// See LICENSE_CELLO file for license and copyright information

/// @file     control_output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @ingroup  Control
/// @brief    Functions implementing CHARM++ output-related functions

// #define DEBUG_OUTPUT

#ifdef DEBUG_OUTPUT
#  define TRACE_OUTPUT(M) CkPrintf ("%d TRACE_OUTPUT %s\n", CkMyPe(), M); fflush(stdout);
#else
#  define TRACE_OUTPUT(M) /*  */
#endif

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void Block::output_enter_ ()
{
  performance_start_(perf_output);
#ifdef NEW_OUTPUT
  new_output_begin_();
#else /* NEW_OUTPUT */
  output_begin_();
#endif  
  performance_stop_(perf_output);
}

//======================================================================
// NEW OUTPUT
//======================================================================

void Block::new_output_begin_ ()
{
  TRACE_OUTPUT("Block::new_output_begin_()");
  cello::simulation() -> new_output_start();
}

//----------------------------------------------------------------------

void Block::new_output_write_block()
{
  TRACE_OUTPUT("Block::new_output_write_block_()");
}

//----------------------------------------------------------------------

void Simulation::new_output_start()
{
  TRACE_OUTPUT("Simulation::new_output_start()");
  if (sync_new_output_start_.next()) {
    TRACE_OUTPUT("Simulation::new_output_start() continuing\n");
    problem()->output_reset();
    problem()->new_output_next(this);
  }
}

//----------------------------------------------------------------------

void Problem::new_output_next(Simulation * simulation) throw()
{
  TRACE_OUTPUT("Problem::new_output_next()");
  int cycle   = simulation->cycle();
  double time = simulation->time();

  Output * output = NULL;

  // Find next schedule output (index_output_ initialized to -1)

  do {

    output = this->output(++index_output_);

  } while (output && ! output->is_scheduled(cycle, time));

  // assert (! output) || ( output->is_scheduled() )

  // print blocks on this process to check block_vec_ array

  Hierarchy * hierarchy = simulation->hierarchy();
  for (int i=0; i<hierarchy->num_blocks(); i++) {
    CkPrintf ("%d Block %d = %s\n",CkMyPe(),i,
	      hierarchy->block(i)->name().c_str());
  }
  if (output != NULL) {

    if (output->is_writer()) {
      TRACE_OUTPUT("Problem::new_output_next: is_writer==true");
    } else {
      TRACE_OUTPUT("Problem::new_output_next: is_writer==false: barrier");
    }

    
    // open file if writer
    output->init();
    output->open();
    output->write_simulation(simulation);
    output->next();


  } else {

    // ...otherwise exit output phase

    TRACE_OUTPUT("Problem::new_output_next(): calling output_exit()");
    simulation->output_exit();

  }
}

//======================================================================
// OLD OUTPUT
//======================================================================

void Block::output_begin_ ()
{
  TRACE_OUTPUT("output_begin_()");
  cello::simulation() -> output_enter();
}

//----------------------------------------------------------------------

void Simulation::output_enter ()
{
  TRACE_OUTPUT("Simulation::output_enter()");

  // Switch from Block to Simulation parallelism
  if (sync_output_begin_.next()) {
    performance_->start_region(perf_output);
    set_phase(phase_output);

    problem()->output_reset();
    problem()->output_next(this);

    performance_->stop_region(perf_output);
  }
}

//----------------------------------------------------------------------

void Problem::output_next(Simulation * simulation) throw()
{
  TRACE_OUTPUT("Problem::output_next()");

  simulation->set_phase(phase_output);

  int cycle   = simulation->cycle();
  double time = simulation->time();

  Output * output;

  // Find next schedule output (index_output_ initialized to -1)

  do {

    output = this->output(++index_output_);

  } while (output && ! output->is_scheduled(cycle, time));

  // assert (! output) || ( output->is_scheduled() )
  
  if (output != NULL) {

    output->next();  // update Output's schedule
    
    // Perform output if any...

    const int stride = output->stride_wait();

    if ((CkMyPe() % stride) == 0) {

      simulation->output_start (index_output_);

    }

  } else {

    // ...otherwise exit output phase

    simulation->output_exit();

  }
}

//----------------------------------------------------------------------

void Simulation::output_start(int index_output)
{
  TRACE_OUTPUT("Simulation::output_start()");
  Output * output = problem()->output(index_output);
  output->init();
  output->open();
  index_output_ = index_output;
  contribute(CkCallback (CkIndex_Simulation::r_output_barrier(NULL),thisProxy));
}

//----------------------------------------------------------------------

void Simulation::r_output_barrier(CkReductionMsg * msg)
{
  delete msg;
  Output * output = problem()->output(index_output_);
  output->write_simulation(this);

  //  ERROR if this is moved here from Output::write_hierarchy()
  //        in OutputCheckpoint since that overrides write_hierarchy()
  //        to not call Output::write_hierarchy_(), which calls
  //        Block::p_output_write()
  
  //  if (CkMyPe() == 0) {
  //    hierarchy_->block_array().p_output_write(index_output,0);
  //  }

}

//----------------------------------------------------------------------

void Block::p_output_write (int index_output, int step)
{
  TRACE_OUTPUT("Simulation::p_output_write()");
  performance_start_ (perf_output);

  FieldDescr    * field_descr    = cello::field_descr();
  ParticleDescr * particle_descr = cello::particle_descr();
  Simulation    * simulation     = cello::simulation();
  Output        * output         = cello::output(index_output);
  Config        * config         = (Config *) cello::config();

  // update derived fields (if any)
  this->compute_derived(config->output_field_list[index_output]);

  output->write_block(this);

  simulation->write_();
  performance_stop_ (perf_output);
}

//----------------------------------------------------------------------

void Simulation::write_()
{
  TRACE_OUTPUT("Simulation::write_()");

  if (sync_output_write_.next()) {

    r_write(NULL);

  }
}

//----------------------------------------------------------------------

void Simulation::r_write(CkReductionMsg * msg)
{
  performance_->start_region(perf_output);
  TRACE_OUTPUT("Simulation::r_write()");
  delete msg;
  problem()->output_wait(this);
  performance_->stop_region(perf_output);
}

//----------------------------------------------------------------------

void Simulation::r_write_checkpoint()
{
  performance_->start_region(perf_output);
  TRACE_OUTPUT("Simulation::r_write_checkpoint()");
  create_checkpoint_link();
  problem()->output_wait(this);
  performance_->stop_region(perf_output);
}

//----------------------------------------------------------------------

void Problem::output_wait(Simulation * simulation) throw()
{
  TRACE_OUTPUT("Problem::output_wait()");
  
  Output * output = this->output(index_output_);

  const int ip = CkMyPe();
  const int np = CkNumPes();
  const int ip_write = output->process_writer();

  if (ip == ip_write) {

    output_write(simulation,0,0);

  } else {

    int n=0;  char * buffer = 0;

    // Copy / alias buffer array of data to send
    output->prepare_remote(&n,&buffer);

    // Send data to writing process
    proxy_simulation[ip_write].p_output_write (n, buffer);

    // Deallocate buffer
    output->cleanup_remote(&n,&buffer);

    output->close();
    output->finalize();
    output_next(simulation);
    
    const int stride = output->stride_wait();
    const int ip_next = ip+1;
    if (ip_next%stride != 0 && ip_next < np) {
      proxy_simulation[ip_next].p_output_start(index_output_);
    }
  }
}

//----------------------------------------------------------------------

void Simulation::p_output_write (int n, char * buffer)
{
  TRACE_OUTPUT("Simulation::p_output_write()");
  problem()->output_write(this,n,buffer); 
}

//----------------------------------------------------------------------

void Problem::output_write 
(
 Simulation * simulation,
 int n, char * buffer
) throw()
{
  TRACE_OUTPUT("Problem::output_write()");

  Output * output = this->output(index_output_);

  if (n != 0) {
    output->update_remote(n, buffer);
  }

  if (output->sync_write()->next()) {

    TRACE_OUTPUT("Problem::output_write(): sync_write()->next() = true");

    output->close();
    output->finalize();
    output_next(simulation);

    const int stride = output->stride_wait();
    const int ip = CkMyPe();
    const int ip_next = ip+1;
    const int np = CkNumPes();
    if (ip_next%stride != 0 && ip_next < np) {
      proxy_simulation[ip_next].p_output_start(index_output_);
    }
  } else {
    TRACE_OUTPUT("Problem::output_write(): sync_write()->next() = false");
  }

}

//----------------------------------------------------------------------

void Simulation::output_exit()
{
  TRACE_OUTPUT("Simulation::output_exit()");

  debug_close();
  debug_open();

  if (CkMyPe() == 0) hierarchy()->block_array().p_output_end();
}

//----------------------------------------------------------------------

void Block::p_output_end()
{
  performance_start_(perf_output);
  TRACE_OUTPUT("Block::p_output_end()");
  performance_stop_(perf_output);
  output_exit_();
}

//======================================================================


