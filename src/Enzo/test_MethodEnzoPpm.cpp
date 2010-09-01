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
#include "user.hpp"
#include "error.hpp"
#include "monitor.hpp"
#include "parameters.hpp"
#include "data.hpp"
#include "parallel.hpp"
#include "test.hpp"
#include "enzo.hpp"

#include "cello_hydro.h"

#include "parallel.def"
#include PARALLEL_CHARM_INCLUDE(test_MethodEnzoPpm.decl.h)

void output_fields(FieldBlock * field_block,
		   int write_count,
		   int field_count,
		   int field_index[],
		   Monitor * monitor)
{

  double map1[] = {1,1,1, 0,0,0};
  FieldDescr * field_descr = field_block->field_descr();
  int nx,ny,nz;
  int gx,gy,gz;
  int mx,my,mz;
  field_block->enforce_boundary(boundary_reflecting);
  field_block->dimensions(&nx,&ny,&nz);
  for (int i = 0; i < field_count; i++) {
    field_descr->ghosts(i,&gx,&gy,&gz);
    mx=nx+2*gx;
    my=ny+2*gy;
    mz=nz+2*gz;
    int index = field_index[i];
    char filename[80];
    std::string field_name = field_descr->field_name(index);
    Scalar * field_values = (Scalar *)field_block->field_values(index);
    sprintf (filename,"ppm-%s-%05d.png",field_name.c_str(),write_count);
    monitor->image (filename, field_values, mx,my,mz, 2, reduce_sum, 0.0,1.0, map1,2);
  }
  
}

PARALLEL_MAIN_BEGIN

{

  // Initialize parallelism

  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  // Create a struct of enzo data (won't work as global data for CHARM++, threading, etc)

  Global * global = new Global;

  unit_init(parallel->rank(), parallel->size());

  // Set necessary parameters for MethodEnzoPpm

  Monitor    * monitor    = global->monitor();
  Parameters * parameters = global->parameters();

  parameters->set_current_group("Physics");
  parameters->set_integer ("dimensions",2);
  parameters->set_scalar  ("gamma",1.4);

  int nx,ny,nz;
  nx=100;
  ny=100;
  nz=1;

  int gx,gy,gz;
  gx=3;
  gy=3;
  gz=0;

  FILE * fp = fopen ("input/test_MethodEnzoPpm.in","r");
  if (fp) {
    parameters->read(fp); // MEMORY LEAK
    fclose(fp);
  } else {
    ERROR_MESSAGE("main","test_MethodEnzoPpm.in file does not exist");
  }

  // Create data and field descriptors and blocks

  DataDescr * data_descr = new DataDescr;
  DataBlock * data_block = new DataBlock;

  FieldDescr * field_descr = data_descr->field_descr();
  FieldBlock * field_block = data_block->field_block();

  // Insert required fields

  int index_density         = field_descr->insert_field("density");
  int index_total_energy    = field_descr->insert_field("total_energy");
  int index_velocity_x      = field_descr->insert_field("velocity_x");
  int index_velocity_y      = field_descr->insert_field("velocity_y");
  int index_internal_energy = field_descr->insert_field("internal_energy");

  EnzoUserDescr    * user_descr = new EnzoUserDescr(global);
  
  UserMethod   * user_method   = user_descr->add_user_method("ppm");
  UserControl  * user_control  = user_descr->set_user_control("ignored");
  UserTimestep * user_timestep = user_descr->set_user_timestep("ignored");

  // Set missing Enzo parameters

  EnzoDescr * enzo = user_descr->enzo();

  enzo->BoundaryRank = 2;
  enzo->BoundaryDimension[0] = nx + 2*gx;
  enzo->BoundaryDimension[1] = ny + 2*gy;
  enzo->BoundaryDimension[2] = nz + 2*gz;

  
  unit_class ("MethodEnzoPpm");

  unit_func("initialize");

  user_control ->initialize(data_descr);
  user_method  ->initialize(data_descr);
  user_timestep->initialize(data_descr);

  // Initialize field_block

  field_block->set_field_descr(field_descr);
  field_block->set_dimensions(nx,ny);
  field_block->set_box_extent(0.0,0.3,0.0,0.3);

  field_descr->set_ghosts (0,gx,gy,gz);
  field_descr->set_ghosts (1,gx,gy,gz);
  field_descr->set_ghosts (2,gx,gy,gz);
  field_descr->set_ghosts (3,gx,gy,gz);
  field_descr->set_ghosts (4,gx,gy,gz);
  
  field_block->allocate_array();
  field_block->allocate_ghosts();
  field_block->clear();

  Scalar *  d = (Scalar * ) field_block->field_values(index_density);
  Scalar * vx = (Scalar * ) field_block->field_values(index_velocity_x);
  Scalar * vy = (Scalar * ) field_block->field_values(index_velocity_y);
  Scalar * te = (Scalar * ) field_block->field_values(index_total_energy);
  //  Scalar * ie = (Scalar * ) field_block->field_values(index_internal_energy);

  double hx = 0.3 / nx;
  double hy = 0.3 / ny;

  int mx = nx + 2*gx;
  int my = ny + 2*gy;

  for (int iy=gy; iy<ny+gy; iy++) {
    double y = (iy - gy + 0.5)*hy;
    for (int ix=gy; ix<nx+gx; ix++) {
      double x = (ix - gx + 0.5)*hx;
      int i = INDEX(ix,iy,0,mx,my);
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

  printf ("Initial\n");

  int indices[] = {index_density,index_velocity_x,index_velocity_y,index_total_energy};

  // Set stopping criteria

  parameters->set_current_group ("Stopping");
  double time_final = parameters->value_scalar("time",2.5);
  int   cycle_final = parameters->value_integer("cycle",1000);
  int  cycle_output = 10;

  double time = 0;
  int cycle = 0;
  //  initialize_implosion(nx+gx);
  user_method  ->initialize_block(data_block);

  // Initial data dump
  output_fields(field_block,cycle,4,indices,monitor);

  while (time < time_final && cycle <= cycle_final) {

    user_control ->initialize_block(data_block);
    user_timestep->initialize_block(data_block);

    field_block->enforce_boundary(boundary_reflecting);

    double dt = user_timestep->compute_block(data_block);

    printf ("cycle = %d  sim-time = %10g dt = %10g\n",
	    cycle,time,dt );

    user_method->advance_block(data_block,time,dt);

    user_timestep->finalize_block(data_block);
    user_control ->finalize_block(data_block);
    user_method  ->finalize_block(data_block);

    ++cycle;
    time += MIN(dt,time_final-time);

    // data dump
    if (cycle % cycle_output == 0) {
      output_fields(field_block,cycle,4,indices,monitor);
    }
    user_timestep->finalize(data_descr);
    user_control ->finalize(data_descr);
  }
  user_method  ->finalize(data_descr);

  unit_finalize();

  delete user_descr;
  delete data_descr;
  delete data_block;
  delete global;
  delete parallel;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(test_MethodEnzoPpm.def.h)

