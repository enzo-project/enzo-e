// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_SimulationMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of SimulationMpi for MPI simulations

#ifndef CONFIG_USE_CHARM

#include "cello.hpp"

#include "simulation.hpp"

//----------------------------------------------------------------------

SimulationMpi::SimulationMpi
(
 const char * parameter_file,
 const GroupProcess * group_process ) throw ()
  : Simulation(parameter_file,group_process)
{
}

//----------------------------------------------------------------------

SimulationMpi::~SimulationMpi() throw()
{
}

//----------------------------------------------------------------------

void SimulationMpi::initialize() throw()
{
  Simulation::initialize();
}

//----------------------------------------------------------------------

void SimulationMpi::run() throw()
{
  
  DEBUG("SimulationMpi::run()");

#ifdef CONFIG_USE_MPI
  ReduceMpi    reduce(group_process_);
#else
  ReduceSerial reduce(group_process_);
#endif

  Problem * problem = Simulation::problem();

  //--------------------------------------------------
  // INITIALIZE FIELDS
  //--------------------------------------------------

  DEBUG("SimulationMpi::run() Initial");

  Initial * initial = problem->initial();

  ItPatch it_patch(hierarchy_);
  Patch * patch;

  while ((patch = ++it_patch)) {

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      initial->enforce_block(block, field_descr_, hierarchy_);

      block->set_cycle(cycle_);
      block->set_time(time_);

      block->initialize();

    }
  }

  // delete parameters in favor of config 
  delete parameters_;
  parameters_ = 0;

  //--------------------------------------------------
  // REFRESH GHOST ZONES AND ENFORCE BOUNDARY CONDITIONS
  //--------------------------------------------------

  DEBUG("SimulationMpi::run() refresh, Boundary ");

  double lower[3];
  hierarchy_->lower(&lower[0], &lower[1], &lower[2]);
  double upper[3];
  hierarchy_->upper(&upper[0], &upper[1], &upper[2]);

  while ((patch = ++it_patch)) {

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      bool is_boundary [3][2];

      block->is_on_boundary (lower,upper,is_boundary);

      refresh_ghost_(block,patch,is_boundary);

      update_boundary_(block,is_boundary);

    }
  }

  //--------------------------------------------------
  // PERFORM SCHEDULED OUTPUT
  //--------------------------------------------------


  DEBUG("SimulationMpi::run() Output");

  scheduled_output();
  
  //--------------------------------------------------
  // EVALUATE INITIAL STOPPING CRITERIA
  //--------------------------------------------------

  Stopping * stopping = problem->stopping();

  int stop_hierarchy = true;

  DEBUG("SimulationMpi::run() Stopping");

  while ((patch = ++it_patch)) {

    int    stop_patch  = true;

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      int cycle = block->cycle();
      int time  = block->time();

      int stop_block = stopping->complete(cycle,time);

      stop_patch = stop_patch && stop_block;
    }
    stop_hierarchy = stop_hierarchy && stop_patch;
  }

  //======================================================================
  // BEGIN MAIN LOOP
  //======================================================================

  while (! stop_hierarchy) {

    //--------------------------------------------------
    // Determine timestep
    //--------------------------------------------------


    Timestep * timestep = problem->timestep();
    Stopping * stopping = problem->stopping();

    double dt_hierarchy = std::numeric_limits<double>::max();

    // Accumulate Patch-local timesteps

    while ((patch = ++it_patch)) {

      double dt_patch = std::numeric_limits<double>::max();

      ItBlock it_block(patch);
      Block * block;

      // Accumulate Block-local timesteps

      while ((block = ++it_block)) {

	double dt_block   = timestep->evaluate(field_descr_,block);
	double time_block = block->time();

	// Reduce timestep to coincide with scheduled output

	int index_output=0;
	while (Output * output = problem->output(index_output++)) {
	  dt_block = output->update_timestep(time_block,dt_block);
	}

	// Reduce timestep to coincide with end of simulation

	double time_stop = stopping->stop_time();
	dt_block = MIN (dt_block, (time_stop - time_block));

	dt_patch = MIN(dt_patch,dt_block);
      }
      dt_hierarchy = MIN(dt_hierarchy, dt_patch);
    }

    dt_hierarchy = reduce.reduce_double(dt_hierarchy,reduce_op_min);

    // Set Simulation timestep

    dt_ = dt_hierarchy;

    // Set Block timesteps

    while ((patch = ++it_patch)) {

      ItBlock it_block(patch);
      Block * block;

      // Accumulate Block-local timesteps

      while ((block = ++it_block)) {
	block->set_dt    (dt_);
      }
    }

    ASSERT("Simulation::run", "dt == 0", dt_ != 0.0);

    // Perform any scheduled output

    scheduled_output();

    monitor_output();

    //--------------------------------------------------
    // REFRESH GHOST ZONES AND ENFORCE BOUNDARY CONDITIONS
    //--------------------------------------------------

    while ((patch = ++it_patch)) {

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	bool is_boundary [3][2];

	block->is_on_boundary (lower,upper,is_boundary);

	refresh_ghost_(block,patch,is_boundary);

	update_boundary_(block,is_boundary);

      }
    }

    //--------------------------------------------------
    // Compute
    //--------------------------------------------------

    while ((patch = ++it_patch)) {

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	int index_method = 0;
	while (Method * method = problem->method(index_method++)) {
	  method -> compute_block (field_descr_,block);
	}

	block->set_cycle (cycle_ + 1);
	block->set_time  (time_ + dt_);

      }
    }
    cycle_ ++;
    time_ += dt_hierarchy;

    //--------------------------------------------------
    // Evaluate stopping criteria
    //--------------------------------------------------

    stop_hierarchy = true;
    while ((patch = ++it_patch)) {
      int stop_patch = true;
      ItBlock it_block(patch);
      Block * block;
      while ((block = ++it_block)) {
	int stop_block = stopping->complete(block->cycle(),block->time());
	stop_patch = stop_patch && stop_block;
      }
      stop_hierarchy = stop_hierarchy && stop_patch;
    }
    stop_hierarchy = reduce.reduce_int(stop_hierarchy,reduce_op_land);


  } // while (! stop_hierarchy)

  //======================================================================
  // END MAIN LOOP
  //======================================================================

  performance_write();
  scheduled_output();
  monitor_output();

}

//----------------------------------------------------------------------

void SimulationMpi::scheduled_output()
{
  size_t index_output = 0;
  while (Output * output = problem()->output(index_output++)) {

    if (output->is_scheduled(cycle_,time_)) {

      output->init();

      output->open();

      output->write_simulation(this);

      //--------------------------------------------------
      int ip       = group_process_->rank();
      int ip_writer = output->process_writer();

      int n=1;  char * buffer = 0;

      if (ip == ip_writer) { // process is writer

	int ip1 = ip+1;
	int ip2 = ip_writer+output->process_stride();

	for (int ip_remote=ip1; ip_remote<ip2; ip_remote++) {

	  // receive size

	  void * handle_recv;
	  handle_recv = group_process_->recv_begin(ip_remote,&n,sizeof(n));
	  group_process_->wait(handle_recv);
	  group_process_->send_end(handle_recv);

	  // allocate buffer

	  buffer = new char [n];

	  // receive buffer

	  handle_recv = group_process_->recv_begin(ip_remote,buffer,n);
	  group_process_->wait(handle_recv);
	  group_process_->recv_end(handle_recv);
	  
	  // update

	  output->update_remote(n,buffer);

	  // deallocate
	  output->cleanup_remote(&n,&buffer);
	}

      } else { // process is not writer

	// send data to writer

	output->prepare_remote(&n,&buffer);

	// send size

	void * handle_send;

	handle_send = group_process_->send_begin(ip_writer,&n,sizeof(n));
	group_process_->wait(handle_send);
	group_process_->send_end(handle_send);

	// send buffer

	handle_send = group_process_->send_begin(ip_writer,buffer,n);
	group_process_->wait(handle_send);
	group_process_->send_end(handle_send);

      }
      //--------------------------------------------------

      output->close();
      output->finalize();
    }
  }
}

//======================================================================

void SimulationMpi::update_boundary_ 
(
 Block * block, 
 bool is_boundary[3][2]
 ) throw()
{
  // Update boundary conditions

  if (dimension_ >= 1) {
    if (is_boundary[axis_x][face_lower]) 
      problem()->boundary()->enforce(field_descr_,block,face_lower,axis_x);
    if (is_boundary[axis_x][face_upper]) 
      problem()->boundary()->enforce(field_descr_,block,face_upper,axis_x);
  }
  if (dimension_ >= 2) {
    if (is_boundary[axis_y][face_lower]) 
      problem()->boundary()->enforce(field_descr_,block,face_lower,axis_y);
    if (is_boundary[axis_y][face_upper]) 
      problem()->boundary()->enforce(field_descr_,block,face_upper,axis_y);
  }
  if (dimension_ >= 3) {
    if (is_boundary[axis_z][face_lower]) 
      problem()->boundary()->enforce(field_descr_,block,face_lower,axis_z);
    if (is_boundary[axis_z][face_upper]) 
      problem()->boundary()->enforce(field_descr_,block,face_upper,axis_z);
  }
}

//----------------------------------------------------------------------

void SimulationMpi::refresh_ghost_ 
(
 Block * block, 
 Patch * patch, 
 bool    is_boundary[3][2]
 ) throw()
{
  // Refresh ghost zones

  int ibx,iby,ibz;
  block->index_patch(&ibx,&iby,&ibz);

  bool periodic = problem()->boundary()->is_periodic();

  // FOLLOWING CODE IS SIMILAR TO Block::x_refresh()

  int nx,ny,nz;
  block->field_block()->size (&nx,&ny,&nz);

  // Determine axes that may be neighbors

  bool axm = (nx > 1) && (periodic || ! is_boundary[axis_x][face_lower]);
  bool axp = (nx > 1) && (periodic || ! is_boundary[axis_x][face_upper]); 
  bool aym = (ny > 1) && (periodic || ! is_boundary[axis_y][face_lower]);
  bool ayp = (ny > 1) && (periodic || ! is_boundary[axis_y][face_upper]);
  bool azm = (nz > 1) && (periodic || ! is_boundary[axis_z][face_lower]);
  bool azp = (nz > 1) && (periodic || ! is_boundary[axis_z][face_upper]);

  int px = (ibx % 2 == 0) ? -1 : 1;
  int py = (iby % 2 == 0) ? -1 : 1;
  int pz = (ibz % 2 == 0) ? -1 : 1;

  if (px == 1) SWAP(axm,axp);
  if (py == 1) SWAP(aym,ayp);
  if (pz == 1) SWAP(azm,azp);

  if (field_descr_->refresh_face(2)) {
    // TRACE3("p %d %d %d",px,py,pz);
    // TRACE6("a %d %d  %d %d  %d %d",axm,axp,aym,ayp,azm,azp);
    if (axm) block->refresh_ghosts(field_descr_,patch,+px,0,0);
    if (axp) block->refresh_ghosts(field_descr_,patch,-px,0,0);
    if (aym) block->refresh_ghosts(field_descr_,patch,0,+py,0);
    if (ayp) block->refresh_ghosts(field_descr_,patch,0,-py,0);
    if (azm) block->refresh_ghosts(field_descr_,patch,0,0,+pz);
    if (azp) block->refresh_ghosts(field_descr_,patch,0,0,-pz);
  }

}

#endif /* ! CONFIG_USE_CHARM */
