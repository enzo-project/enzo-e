// $Id: test_enzo_method_ppml.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     test_enzo_method_ppml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Apr  1 16:19:18 PDT 2010
/// @brief    Unit tests for the EnzoMethodPpml class

#include <stdio.h>
#include <mpi.h>

#include "parallel.hpp"
#include "user.hpp"
#include "test.hpp"

int main (int argc, char ** argv)
{

  // Initialize parallelism

  Parallel * parallel = Parallel::instance();
  parallel->initialize(&argc,&argv);

  unit_init(parallel->process_rank(), parallel->process_count());

  unit_class ("MethodEnzoPpml");
  MethodEnzoPpml ppml;

  unit_func("initialize");
  ppml.initialize();
  unit_assert(true);

  unit_func("advance_block");
  ppml.advance_block();
  unit_assert(false);

  unit_func("refresh_face");
  ppml.refresh_face();
  unit_assert(false);

  unit_finalize();

}
