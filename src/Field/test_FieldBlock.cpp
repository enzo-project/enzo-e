// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Block.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the FieldBlock class

#include "cello.h"

#include "error.hpp"
#include "test.hpp"
#include "field.hpp"

int main()
{

  //----------------------------------------------------------------------
  unit_init();
  //----------------------------------------------------------------------

  FieldDescr field_descr;

  int index_density      = field_descr.insert_field("density");
  int index_velocity_x   = field_descr.insert_field("velocity_x");
  int index_velocity_y   = field_descr.insert_field("velocity_y");
  int index_velocity_z   = field_descr.insert_field("velocity_z");
  int index_total_energy = field_descr.insert_field("total_energy");

  field_descr.set_precision(index_density,     precision_single);
  field_descr.set_precision(index_velocity_x,  precision_double);
  field_descr.set_precision(index_velocity_y,  precision_quadruple);
  field_descr.set_precision(index_velocity_z,  precision_double);
  field_descr.set_precision(index_total_energy,precision_half);

  field_descr.set_ghosts(index_density,      1,1,1);
  field_descr.set_ghosts(index_velocity_x,   2,2,2);
  field_descr.set_ghosts(index_velocity_y,   3,2,1);
  field_descr.set_ghosts(index_velocity_z,   1,2,3);
  field_descr.set_ghosts(index_total_energy, 5,5,5);

  field_descr.set_centering(index_velocity_x, false, true, true);
  field_descr.set_centering(index_velocity_y, true, false, true);
  field_descr.set_centering(index_velocity_z, true, true,  false);

  //  printf ("sizeof(half) = %d\n",sizeof(float16));
  printf ("sizeof(single) = %d\n",sizeof(float));
  printf ("sizeof(double) = %d\n",sizeof(double));
  printf ("sizeof(extended) = %d\n",sizeof(long double));

  //----------------------------------------------------------------------
  unit_class ("FieldBlock");
  //----------------------------------------------------------------------

  FieldBlock field_block;

  //----------------------------------------------------------------------

  unit_func("field_descr");

  field_block.set_field_descr(&field_descr);
  unit_assert (field_block.field_descr() == & field_descr);

  //----------------------------------------------------------------------

  unit_func("dimensions");

  field_block.set_dimensions(4,5,6);
  int dimensions[3];
  field_block.dimensions(&dimensions[0],&dimensions[1],&dimensions[2]);
  unit_assert(dimensions[0]==4 && dimensions[1]==5 && dimensions[2]==6);

  field_block.set_dimensions(5,3,4);
  field_block.dimensions(&dimensions[0],&dimensions[1],&dimensions[2]);
  unit_assert(dimensions[0]==5 && dimensions[1]==3 && dimensions[2]==4);

  //----------------------------------------------------------------------
  // allocate / deallocate
  //----------------------------------------------------------------------

  unit_func("array_allocated");
  unit_assert( ! field_block.array_allocated());

  // Allocate

  unit_func("allocate_array");
  field_block.allocate_array();
  unit_assert(field_block.array() != 0);

  unit_func("array_allocated");
  unit_assert( field_block.array_allocated());

  // Deallocate

  unit_func("deallocate_array");
  field_block.deallocate_array();
  unit_assert(field_block.array() == 0);

  unit_func("array_allocate");
  unit_assert( ! field_block.array_allocated());

  // Allocate

  unit_func("allocate_array");
  field_block.allocate_array();
  unit_assert(field_block.array() != 0);
  
  unit_func("array_allocated");
  unit_assert( field_block.array_allocated());

  
  //----------------------------------------------------------------------

  unit_func("field_values");
  
  float * values_density = (float *) field_block.field_values(index_density);
  

  unit_assert(false);
  unit_func("field_unknowns");
  unit_assert(false);

  //----------------------------------------------------------------------
  unit_func("box_extent");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("cell_width");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_func("clear");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_func("ghosts_allocated");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("allocate_ghosts");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("deallocate_ghosts");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_func("split");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("merge");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_func("read");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("write");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_finalize();
  //----------------------------------------------------------------------
}
