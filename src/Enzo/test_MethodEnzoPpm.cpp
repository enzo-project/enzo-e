// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_MethodEnzoPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:04:03 PST 2008
/// @brief    Program implementing unit tests for the MethodEnzoPpm class

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cello.hpp"
#include "cello_hydro.h"

#include "method.hpp"
#include "error.hpp"
#include "monitor.hpp"
#include "parameters.hpp"
#include "parallel.hpp"
#include "test.hpp"
#include "enzo.hpp"
#include "simulation.hpp"
#include "parallel.def"

#include PARALLEL_CHARM_INCLUDE(test_MethodEnzoPpm.decl.h)

//----------------------------------------------------------------------

void output_images(Monitor * monitor,
		   int cycle,
		   FieldBlock * field_block,
		   int field_count,
		   int field_index[],
		   int cycle_image = 1)
{

  if (! (cycle_image && cycle % cycle_image == 0)) return;

  FieldDescr * field_descr = field_block->field_descr();
  int nx,ny,nz;
  int gx,gy,gz;
  int mx,my,mz;
  field_block->enforce_boundary(boundary_reflecting);
  field_block->size(&nx,&ny,&nz);
  for (int i = 0; i < field_count; i++) {
    field_descr->ghosts(i,&gx,&gy,&gz);
    mx=nx+2*gx;
    my=ny+2*gy;
    mz=nz+2*gz;
    int index = field_index[i];
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
		 int field_count,
		 int field_index[],
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
  for (int i = 0; i < field_count; i++) {

    // Get field's ghost zone depth
    int gx,gy,gz;
    field_descr->ghosts(i,&gx,&gy,&gz);

    // Get field's total size including ghosts

    int mx,my,mz;
    mx=nx+2*gx;
    my=ny+2*gy;
    mz=nz+2*gz;

    // Get field's index

    int index = field_index[i];

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

  Global     * global = new Global;

  Monitor    * monitor    = global->monitor();
  Parameters * parameters = global->parameters();
  FileHdf5   hdf5;

  // Read in parameters

  FILE * file_pointer;
  file_pointer = fopen ("input/test_MethodEnzoPpm.in","r");
  parameters->read(file_pointer); // MEMORY LEAK
  fclose(file_pointer);

  // Create top-level Simulation object

  EnzoSimulation simulation(global);

  simulation.initialize();

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // Create and initialize data descriptor

  // Mesh::Mesh()
  parameters->set_current_group ("Mesh");

  // Mesh size

  // Mesh::Mesh()
  int nx = parameters->list_value_integer(0,"root_size",1);
  int ny = parameters->list_value_integer(1,"root_size",1);
  int nz = parameters->list_value_integer(2,"root_size",1);

  // Block size

  // Mesh::Mesh()
  int mx = parameters->list_value_integer(0,"block_size",1);
  int my = parameters->list_value_integer(1,"block_size",1);
  int mz = parameters->list_value_integer(2,"block_size",1);

  // Block count

//   printf ("%d %d %d  %d %d %d\n",nx,ny,nz,mx,my,mz);

  int kx = int (ceil (1.0*nx/mx));
  int ky = int (ceil (1.0*ny/my));
  int kz = int (ceil (1.0*nz/mz));

  DataBlock ** data_block_array = new DataBlock * [kx*ky*kz];
  for (int k=0; k<kx*ky*kz; k++) {
    data_block_array[k] = new DataBlock;
  }

  ASSERT("main","Assuming only one data_block",kx*ky*kz == 1);

  DataBlock * data_block = data_block_array[0];

  FieldBlock * field_block = data_block->field_block();

  // Set missing Enzo parameters

  EnzoDescr * enzo = simulation.enzo();

  // Initialize field_block

  unit_class ("FieldBlock");
  unit_func("initialize");

  DataDescr  * data_descr  = simulation.data_descr();
  FieldDescr * field_descr = data_descr->field_descr();

  parameters->set_current_group ("Domain");
  
  double lower[3];
  double upper[3];

  lower[0] = parameters->list_value_scalar(0,"extent",0.0);
  lower[1] = parameters->list_value_scalar(2,"extent",0.0);
  lower[2] = parameters->list_value_scalar(4,"extent",0.0);

  upper[0] = parameters->list_value_scalar(1,"extent",1.0);
  upper[1] = parameters->list_value_scalar(3,"extent",1.0);
  upper[2] = parameters->list_value_scalar(5,"extent",1.0);

  field_block->set_field_descr(field_descr);
  field_block->set_size(nx,ny);
  field_block->set_extent(0.0,0.3,0.0,0.3);

  int gx=3,gy=3,gz=0;
  field_descr->set_ghosts (0,gx,gy,gz);
  field_descr->set_ghosts (1,gx,gy,gz);
  field_descr->set_ghosts (2,gx,gy,gz);
  field_descr->set_ghosts (3,gx,gy,gz);
  field_descr->set_ghosts (4,gx,gy,gz);
  
  field_block->allocate_array();
  field_block->allocate_ghosts();
  field_block->clear();

  int index_density = 0;
  int index_velocity_x = 1;
  int index_velocity_y = 2;
  int index_total_energy = 3;
  int index_internal_energy = 4;
  
  Scalar *  d = (Scalar * ) field_block->field_values(index_density);
  Scalar * vx = (Scalar * ) field_block->field_values(index_velocity_x);
  Scalar * vy = (Scalar * ) field_block->field_values(index_velocity_y);
  Scalar * te = (Scalar * ) field_block->field_values(index_total_energy);
  //  Scalar * ie = (Scalar * ) field_block->field_values(index_internal_energy);

  double hx = 0.3 / nx;
  double hy = 0.3 / ny;

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;

  for (int iy=gy; iy<ny+gy; iy++) {
    double y = (iy - gy + 0.5)*hy;
    for (int ix=gy; ix<nx+gx; ix++) {
      double x = (ix - gx + 0.5)*hx;
      int i = INDEX(ix,iy,0,ngx,ngy);
      if (x + y < 0.1517) {
	d[i]  = 0.125;
	vx[i] = 0;
	vy[i] = 0;
	te[i] = 0.14 / ((enzo->Gamma - 1.0) * d[i]);
      } else {
	d[i]  = 1.0;
	vx[i] = 0;
	vy[i] = 0;
	te[i] = 1.0 / ((enzo->Gamma - 1.0) * d[i]);
      }
    }
  }

  int indices[] = {index_density,
		   index_velocity_x,
		   index_velocity_y,
		   index_total_energy};

  // Initialize stopping criteria

  parameters->set_current_group ("Stopping");

  double time_final = parameters->value_scalar("time",2.5);
  int   cycle_final = parameters->value_integer("cycle",1000);

  double time = 0;
  double dt = 0;
  int cycle = 0;

  // Initialize Output parameters

  parameters->set_current_group ("Output");

  int  cycle_dump    = parameters->value_integer("cycle_dump",10);

  // Initialize monitoring parameters

  parameters->set_current_group ("Monitor");

  int  cycle_image    = parameters->value_integer("cycle_image",10);
  int  cycle_progress = parameters->value_integer("cycle_progress",1);

  // Initial progress and image monitoring

  output_progress(monitor,cycle,time,dt);
  output_images(monitor,cycle,field_block,4,indices);

  Papi papi;
  Timer timer;
  
  timer.start();
  papi.start();
  while (time < time_final && cycle <= cycle_final) {

    simulation.initialize_block(data_block);

    field_block->enforce_boundary(boundary_reflecting);

    dt = simulation.method_timestep()->compute_block(data_block);

    output_images  (monitor,cycle,field_block,4,indices,cycle_image);
    output_progress(monitor,cycle,time,dt,cycle_progress);
    output_dump    (hdf5,cycle,field_block,4,indices,cycle_dump);

    simulation.method_hyperbolic(0)->advance_block(data_block,time,dt);

    simulation.finalize_block(data_block);

    ++cycle;
    time += MIN(dt,time_final-time);

  }

  papi.stop();
  timer.stop();

  output_images  (monitor,cycle_final,field_block,4,indices);
  output_progress(monitor,cycle_final,time,dt);

  printf ("Time real = %f\n",papi.time_real());
  printf ("Time proc = %f\n",papi.time_proc());
  printf ("Flop count = %lld\n",papi.flop_count());
  printf ("GFlop rate = %f\n",papi.flop_rate()*1e-9);

  printf ("Timer time = %f\n",timer.value());

  simulation.finalize();

  unit_finalize();

  delete data_block;
  delete global;
  delete parallel;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(test_MethodEnzoPpm.def.h)

