
// $Id: test_MethodEnzoPpm.cpp 1552 2010-06-09 02:50:30Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_MethodEnzoPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:04:03 PST 2008
/// @brief    Program implementing unit tests for the MethodEnzoPpm class

#include <mpi.h>
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

int main(int argc, char **argv)
{

  // Initialize parallelism

  Parallel * parallel = Parallel::instance();
  parallel->initialize(&argc,&argv);

  Parameters * parameters = Parameters::instance();

  unit_init(parallel->process_rank(), parallel->process_count());

  // Create data and field descriptors and blocks

  FieldDescr * field_descr = new FieldDescr;
  FieldBlock * field_block = new FieldBlock;

  DataDescr * data_descr = new DataDescr (field_descr);
  DataBlock * data_block = new DataBlock (field_block);

  // Insert required fields

  int index_density         = field_descr->insert_field("density");
  int index_total_energy    = field_descr->insert_field("total_energy");
  int index_internal_energy = field_descr->insert_field("internal_energy");
  int index_velocity_x      = field_descr->insert_field("velocity_x");
  int index_velocity_y      = field_descr->insert_field("velocity_y");

  // Initialize field_block

  int nx,ny,nz;
  nx=100;
  ny=100;
  nz=1;
  field_block->set_field_descr(field_descr);
  field_block->set_dimensions(nx,ny);
  field_block->set_box_extent(0.0,0.3,0.0,0.3);

  int gx,gy,gz;
  gx=3;
  gy=3;
  gz=0;
  field_descr->set_ghosts (0,gx,gy,gz);
  field_descr->set_ghosts (1,gx,gy,gz);
  field_descr->set_ghosts (2,gx,gy,gz);
  field_descr->set_ghosts (3,gx,gy,gz);
  field_descr->set_ghosts (4,gx,gy,gz);
  
  field_block->allocate_array();
  field_block->allocate_ghosts();
  field_block->clear();

  Scalar * d = (Scalar * ) field_block->field_values(index_density);
  Scalar * vx = (Scalar * ) field_block->field_values(index_velocity_x);
  Scalar * vy = (Scalar * ) field_block->field_values(index_velocity_y);
  Scalar * te = (Scalar * ) field_block->field_values(index_internal_energy);

  double hx = 0.3 / nx;
  double hy = 0.3 / ny;

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  for (int iy=gy; iy<ny+gy; iy++) {
    double y = (iy - gy + 0.5)*hy;
    for (int ix=gy; ix<nx+gx; ix++) {
      double x = (ix - gx + 0.5)*hx;
      int i = INDEX(ix,iy,0,mx,my);
      if (x + y < 0.1517) {
	d[i]  = 0.125;
	vx[i] = 0;
	vy[i] = 0;
	te[i] = 0.14 / ((1.4 - 1.0) * d[i]);
      } else {
	d[i]  = 1.0;
	vx[i] = 0;
	vy[i] = 0;
	te[i] = 1.0 / ((1.4 - 1.0) * d[i]);
      }
    }
  }

  // Set necessary parameters for MethodEnzoPpm

  parameters->set_current_group("Physics");
  parameters->set_integer ("dimensions",2);

  unit_class ("MethodEnzoPpm");
  MethodEnzoPpm ppm;

  unit_func("initialize_method");
  ppm.initialize_method(data_descr);
  unit_assert(true);

  Monitor * monitor = Monitor::instance();

  double map1[] = {0,0,0, 1,1,1};

  field_block->enforce_boundary(boundary_reflecting);

  monitor->image ("ppm-density-0.png",
		  (Scalar *)field_block->field_values(index_density),
		  mx,my,mz,
		  0,  0,  0,
		  mx,my,mz,
		  2,
		  reduce_sum,
		  0.0,1.0,
		  map1,2);

  unit_func("initialize_block");
  ppm.initialize_block(data_block);
  unit_assert(true);

  double t = 0;
  double dt = 0.000239579;

  unit_func("advance_block");
  ppm.advance_block(data_block,t,dt);
  unit_assert(false);

  unit_func("finalize_block");
  ppm.finalize_block(data_block);
  unit_assert(false);

  monitor->image ("ppm-density-1.png",
		  (Scalar *)field_block->field_values(index_density),
		  mx,my,mz,
		  0,  0,  0,
		  mx,my,mz,
		  2,
		  reduce_sum,
		  0.0,1.0,
		  map1,2);

  unit_func("refresh_face");
  ppm.refresh_face();
  unit_assert(false);

  unit_func("finalize_method");
  ppm.finalize_method(data_descr);
  unit_assert(false);

  unit_finalize();

  parallel->finalize();
}

