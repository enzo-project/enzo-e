// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file     test_MethodEnzoPpml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Apr  1 16:19:18 PDT 2010
/// @brief    Unit tests for the EnzoMethodPpml class

#include <stdio.h>
#include "data.hpp"
#include "parallel.hpp"
#include "user.hpp"
#include "test.hpp"

int main (int argc, char ** argv)
{

  // Initialize parallelism

  Mpi::init(&argc,&argv);

  Parallel * parallel = Parallel::instance();

  unit_init(parallel->process_rank(), parallel->process_count());

  // Create a struct of enzo data (won't work as global data for CHARM++, threading, etc)

  DataDescr * data_descr = new DataDescr(new FieldDescr);
  DataBlock * data_block = new DataBlock(new FieldBlock);

  unit_class ("MethodEnzoPpml");
  MethodEnzoPpml ppml;

  unit_func("initialize_block");
  ppml.initialize_block(data_block);
  unit_assert(true);

  double t = 0;
  double dt = 0.1;

  unit_func("advance_block");
  ppml.advance_block(data_block,t,dt);
  unit_assert(false);

  unit_func("finalize_block");
  ppml.finalize_block(data_block);
  unit_assert(false);

  unit_func("refresh_face");
  ppml.refresh_block();
  unit_assert(false);

  unit_func("finalize_method");
  ppml.finalize(data_descr);
  unit_assert(false);

  unit_finalize();

  Mpi::finalize();

}
