// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:04:03 PST 2008
/// @todo     Move initialization into Simulation / EnzoSimulation
/// @brief    Program implementing unit tests for the EnzoMethodPpm class

#include "cello.hpp"

#include "enzo.hpp"

#include "test.hpp"

#include PARALLEL_CHARM_INCLUDE(test_EnzoMethodPpm.decl.h)

//----------------------------------------------------------------------

void output_images(Monitor * monitor,
		   int cycle,
		   FieldBlock * field_block,
		   int cycle_image = 1)
{

  if (! (cycle_image && cycle % cycle_image == 0)) return;

  FieldDescr * field_descr = field_block->field_descr();
  int nx,ny,nz;
  int gx,gy,gz;
  int mx,my,mz;
  field_block->enforce_boundary(boundary_reflecting);
  field_block->size(&nx,&ny,&nz);
  int count = field_descr->field_count();
  for (int index = 0; index < count; index++) {
    field_descr->ghosts(index,&gx,&gy,&gz);
    mx=nx+2*gx;
    my=ny+2*gy;
    mz=nz+2*gz;
    char filename[80];
    std::string field_name = field_descr->field_name(index);
    Scalar * field_values = (Scalar *)field_block->field_values(index);
    sprintf (filename,"ppm-%s-%05d.png",field_name.c_str(),cycle);
    monitor->image (filename, field_values, mx,my,mz, 2, reduce_sum, 0.0, 1.0);
  }
}

//----------------------------------------------------------------------

void output_dump(FileHdf5 & hdf5,
		 int cycle,
		 FieldBlock * field_block,
		 int cycle_dump = 1)
{

  // Exit if we don't dump data this cycle
  if (! (cycle_dump && cycle % cycle_dump == 0)) return;

  // Refresh boundary conditions 
  // (should have check to not do it more than once)

  FieldDescr * field_descr = field_block->field_descr();
  field_block->enforce_boundary(boundary_reflecting);

  // Open file

  char filename[80];
  sprintf (filename,"ppm-%05d.h5",cycle);
  hdf5.file_open (filename,"w");

  // Get block size
  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  // Loop over fields in block

  int count = field_descr->field_count();

  for (int index = 0; index < count; index++) {

    // Get field's ghost zone depth
    int gx,gy,gz;
    field_descr->ghosts(index,&gx,&gy,&gz);

    // Get field's total size including ghosts

    int mx,my,mz;
    mx=nx+2*gx;
    my=ny+2*gy;
    mz=nz+2*gz;

    // Get field's name for filename

    std::string field_name = field_descr->field_name(index);

    // Write data_block, including ghosts
    // (should include option to omit ghosts)
    // (what if ghosts aren't allocated in general?)

    // prepare to write the field to the file

    hdf5.dataset_open_write (field_name,
			     field_descr->precision(index),
			     mx,my,mz);

    // write the field to the file

    hdf5.write((char *)field_block->field_values(index),
		field_descr->precision(index));

    // close the field in the file
    hdf5.dataset_close ();
  }

  // close the file

  hdf5.file_close();

}

//----------------------------------------------------------------------

void output_progress (Monitor * monitor,
		      int cycle,
		      double time,
		      double dt,
		      int cycle_progress = 1
		      )
{
  if (! (cycle_progress && cycle % cycle_progress == 0)) return;
  char buffer[100];
  sprintf (buffer," cycle = %05d  sim-time = %10.8f  dt = %10.8f",
	   cycle,time,dt);
  monitor->print (buffer);
}


PARALLEL_MAIN_BEGIN

{

  // Initialize parallelism

  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();
  unit_init(parallel->rank(), parallel->size());

  // Initialize cross-cutting components

  Error     * error    = new Error;
  Monitor   * monitor  = new Monitor;

  // Create top-level Simulation object

  Simulation * simulation = new EnzoSimulation(error,monitor);

  FILE * fp = fopen ("input/test_EnzoMethodPpm.in","r");

  simulation->initialize(fp);

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  //   FileHdf5   hdf5;

  // output_progress(monitor,cycle,time,dt);
  // output_images(monitor,cycle,field_block);

  // Papi papi;
  // Timer timer;
  
  // timer.start();
  // papi.start();

  // while (time < time_final && cycle <= cycle_final) {

  //   simulation->initialize_block(data_block);

  //   field_block->enforce_boundary(boundary_reflecting);

  //   dt = simulation->imestep()->compute_block(data_block);

  //   output_images  (monitor,cycle,field_block,cycle_image);
  //   output_progress(monitor,cycle,time,dt,cycle_progress);
  //   output_dump    (hdf5,cycle,field_block,cycle_dump);

  //   simulation->method(0)->compute_block(data_block,time,dt);

  //   simulation->finalize_block(data_block);

  //   ++cycle;
  //   time += MIN(dt,time_final-time);

  // }

  // papi.stop();
  // timer.stop();

  // output_images  (monitor,cycle_final,field_block);
  // output_progress(monitor,cycle_final,time,dt);

  // printf ("Time real = %f\n",papi.time_real());
  // printf ("Time proc = %f\n",papi.time_proc());
  // printf ("Flop count = %lld\n",papi.flop_count());
  // printf ("GFlop rate = %f\n",papi.flop_rate()*1e-9);

  // printf ("Timer time = %f\n",timer.value());

  simulation->finalize();

  unit_finalize();

  delete error;
  delete monitor;
  delete simulation;

  fclose (fp);

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(test_EnzoMethodPpm.def.h)

