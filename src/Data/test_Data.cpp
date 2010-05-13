// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_data.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the DataBlock class

#include "cello.h"

#include "error.hpp"
#include "test.hpp"
#include "data.hpp"
#include "particle.hpp"
#include "field.hpp"

int main()
{

  unit_class ("DataDescr");

  DataDescr data;

  unit_func("set_dimension");
  data.set_dimension(2);
  unit_assert (data.dimension() == 2);

  unit_func("field_insert");
  unit_assert(data.field_count() == 0);
  int index_density = data.field_insert("density");
  unit_assert(data.field_count() == 1);
  int index_velocity_x = data.field_insert("velocity_x");
  unit_assert(data.field_count() == 2);
  int index_velocity_y = data.field_insert("velocity_y");
  unit_assert(data.field_count() == 3);

  unit_func("field_descr");
  FieldDescr * density    = data.field_descr (index_density);
  FieldDescr * velocity_x = data.field_descr (index_velocity_x);
  FieldDescr * velocity_y = data.field_descr (index_velocity_y);

  unit_assert(density->name()    == "density");
  unit_assert(velocity_x->name() == "velocity_x");
  unit_assert(velocity_y->name() == "velocity_y");

  unit_func ("field_index");
  unit_assert(data.field_index("density")    == index_density);
  unit_assert(data.field_index("velocity_x") == index_velocity_x);
  unit_assert(data.field_index("velocity_y") == index_velocity_y);

  unit_func ("field_name");
  unit_assert(data.field_name(index_density)    == "density");
  unit_assert(data.field_name(index_velocity_x) == "velocity_x");
  unit_assert(data.field_name(index_velocity_y) == "velocity_y");

  unit_func("set_dimension");
  unit_assert(density->   dimension()==2);
  unit_assert(velocity_x->dimension()==2);
  unit_assert(velocity_y->dimension()==2);

}
