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

  field_descr.insert_field("density");
  field_descr.insert_field("velocity_x");
  field_descr.insert_field("velocity_y");
  field_descr.insert_field("velocity_z");
  field_descr.insert_field("total_energy");

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
  unit_func("index_range");
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
