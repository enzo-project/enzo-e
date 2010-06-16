
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

#include "user.hpp"
#include "error.hpp"
#include "data.hpp"
#include "parallel.hpp"
#include "test.hpp"

int main(int argc, char **argv)
{

  // Initialize parallelism

  Parallel * parallel = Parallel::instance();
  parallel->initialize(&argc,&argv);

  unit_init(parallel->process_rank(), parallel->process_count());

  
  DataDescr * data_descr = new DataDescr (new FieldDescr);
  DataBlock * data_block = new DataBlock (new FieldBlock);

  unit_class ("MethodEnzoPpm");
  MethodEnzoPpm ppm;

  unit_func("initialize_method");
  ppm.initialize_method(data_descr);
  unit_assert(true);

  unit_func("initialize_block");
  ppm.initialize_block(data_block);
  unit_assert(true);

  double t = 0;
  double dt = 0.1;

  unit_func("advance_block");
  ppm.advance_block(data_block,t,dt);
  unit_assert(false);

  unit_func("finalize_block");
  ppm.finalize_block(data_block);
  unit_assert(false);

  unit_func("refresh_face");
  ppm.refresh_face();
  unit_assert(false);

  unit_func("finalize_method");
  ppm.finalize_method(data_descr);
  unit_assert(false);

  unit_finalize();

  parallel->finalize();
}

