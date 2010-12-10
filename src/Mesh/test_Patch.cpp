// $Id: test_Mesh.cpp 1791 2010-10-29 22:54:06Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Mesh class
 
#include <stdio.h>
#include <string>

#include "cello.hpp"
#include "test.hpp"
#include "mesh.hpp"

#include "parallel.def"

#include PARALLEL_CHARM_INCLUDE(test_Mesh.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();
  unit_class ("Mesh");
  unit_func("Mesh");


  void set_data_descr (DataDescr * data_descr) throw();


  DataDescr * data_descr () throw();



  void set_patch_size (int npx, int npy, int npz) throw();



  void patch_size (int * npx, int * npy=1, int * npz=1) throw();



  void set_layout (Layout * layout) throw();



  Layout * layout () throw();



  void set_extents (double xm, double xp,
		    double ym, double yp,
		    double zm, double zp) throw();



  void extents (double * xm, double * xp,
		double * ym=0, double * yp=0,
		double * zm=0, double * zp=0) throw();

  
  void p_evolve();


  
void allocate_blocks() ();

/// Deallocate local blocks
void deallocate_blocks() ();


  int block_count()  throw();
  

  DataBlock * block(int i) throw();

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Mesh.def.h)
