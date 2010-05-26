// $Id: test_field.cpp 1412 2010-05-05 00:10:09Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Field.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Test program for the FieldDescr class

#include "cello.h"

#include "test.hpp"
#include "field.hpp"

int main()
{

  unit_init();
  unit_class ("FieldDescr");

  FieldDescr * field_descr = 0;

  unit_func ("FieldDescr");
  field_descr = new FieldDescr;
  unit_assert(field_descr != 0);

  unit_func("~FieldDescr");
  delete field_descr;
  unit_assert(true);

  field_descr = new FieldDescr;

  // Fields

  unit_func("insert_field");
  unit_assert(field_descr->field_count()==0);
  field_descr->insert_field("density");
  unit_assert(field_descr->field_count()==1);
  field_descr->insert_field("velocity_x");
  unit_assert(field_descr->field_count()==2);
  field_descr->insert_field("velocity_y");
  unit_assert(field_descr->field_count()==3);
  field_descr->insert_field("velocity_z");
  unit_assert(field_descr->field_count()==4);
  field_descr->insert_field("total_energy");

  unit_func("field_count");
  unit_assert(field_descr->field_count()==5);

  unit_func("field_id");

  int field_density      = field_descr->field_id("density");
  int field_velocity_x   = field_descr->field_id("velocity_x");
  int field_velocity_y   = field_descr->field_id("velocity_y");
  int field_velocity_z   = field_descr->field_id("velocity_z");
  int field_total_energy = field_descr->field_id("total_energy");

  unit_assert(field_descr->field_id("density")      == field_density);
  unit_assert(field_descr->field_id("velocity_x")   == field_velocity_x);
  unit_assert(field_descr->field_id("velocity_y")   == field_velocity_y);
  unit_assert(field_descr->field_id("velocity_z")   == field_velocity_z);
  unit_assert(field_descr->field_id("total_energy") == field_total_energy);

  unit_func("field_name");

  unit_assert(field_descr->field_name(field_density)      == "density");
  unit_assert(field_descr->field_name(field_velocity_x)   == "velocity_x");
  unit_assert(field_descr->field_name(field_velocity_y)   == "velocity_y");
  unit_assert(field_descr->field_name(field_velocity_z)   == "velocity_z");
  unit_assert(field_descr->field_name(field_total_energy) == "total_energy");

  // Groups

  unit_func("insert_group");
  unit_assert(field_descr->group_count()==0);
  field_descr->insert_group("density");
  unit_assert(field_descr->group_count()==1);
  field_descr->insert_group("vector");
  unit_func("group_count");
  unit_assert(field_descr->group_count()==2);

  unit_func("group_id");

  int group_density = field_descr->group_id("density");
  int group_vector   = field_descr->group_id("vector");

  unit_assert(field_descr->group_id("density") == group_density);
  unit_assert(field_descr->group_id("vector")   == group_vector);

  unit_func("group_name");

  unit_assert(field_descr->group_name(group_density) == "density");
  unit_assert(field_descr->group_name(group_vector)    == "vector");

  // Fields and groups

  unit_func("set_field_in_group");
  field_descr->set_field_in_group(field_density,  group_density);
  field_descr->set_field_in_group(field_velocity_x,group_vector);
  field_descr->set_field_in_group(field_velocity_y,group_vector);
  field_descr->set_field_in_group(field_velocity_z,group_vector);

  unit_assert(field_descr->field_in_group(field_density,   group_density));
  unit_assert(field_descr->field_in_group(field_velocity_x,group_vector));
  unit_assert(field_descr->field_in_group(field_velocity_y,group_vector));
  unit_assert(! field_descr->field_in_group(field_velocity_y,  group_density));
  unit_assert(! field_descr->field_in_group(field_total_energy,group_density));
  unit_assert(! field_descr->field_in_group(field_density,     group_vector));
  unit_assert(! field_descr->field_in_group(field_total_energy,group_vector));


  //----------------------------------------------------------------------
  // Global attributes
  //----------------------------------------------------------------------

  // (set and reset in case test value is a default)

  field_descr->set_alignment(8);
  field_descr->set_padding(64);
  field_descr->set_courant(0.5);
  
  unit_func("alignment");
  unit_assert(field_descr->alignment() == 8);
  unit_func("padding");
  unit_assert(field_descr->padding() == 64);
  unit_func("courant");
  unit_assert(field_descr->courant() == 0.5);

  field_descr->set_alignment(4);
  field_descr->set_padding(32);
  field_descr->set_courant(0.75);
  
  unit_func("alignment");
  unit_assert(field_descr->alignment() == 4);
  unit_func("padding");
  unit_assert(field_descr->padding() == 32);
  unit_func("courant");
  unit_assert(field_descr->courant() == 0.75);
  
  //----------------------------------------------------------------------
  // Field attributes
  //----------------------------------------------------------------------

  // Precision

  unit_func("precision");

  field_descr->set_precision(field_density,    precision_single);
  field_descr->set_precision(field_velocity_x, precision_double);
  field_descr->set_precision(field_velocity_y, precision_double);
  field_descr->set_precision(field_velocity_z, precision_double);

  unit_assert(field_descr->precision(field_density)      == precision_single);
  unit_assert(field_descr->precision(field_velocity_x)   == precision_double);
  unit_assert(field_descr->precision(field_velocity_y)   == precision_double);
  unit_assert(field_descr->precision(field_velocity_z)   == precision_double);
  unit_assert(field_descr->precision(field_total_energy) == precision_default);

  // Centering

  unit_func("centering");

  field_descr->set_centering(field_velocity_x, false, true, true);
  field_descr->set_centering(field_velocity_y, true, false, true);
  field_descr->set_centering(field_velocity_z, true,  true, false);

  bool cx, cy, cz;

  field_descr->centering(field_density, &cx, &cy, &cz);
  unit_assert(cx && cy && cz);

  field_descr->centering(field_velocity_x, &cx, &cy, &cz);
  unit_assert(! cx &&   cy &&   cz);

  field_descr->centering(field_velocity_y, &cx, &cy, &cz);
  unit_assert(  cx && ! cy &&   cz);

  field_descr->centering(field_velocity_z, &cx, &cy, &cz);
  unit_assert(  cx &&   cy && ! cz);
  
  // Ghost zone depth

  unit_func("ghosts");

  field_descr->set_ghosts(field_density, 3, 3, 3);
  field_descr->set_ghosts(field_velocity_x, 1, 0, 0);
  field_descr->set_ghosts(field_velocity_y, 0, 1, 0);
  field_descr->set_ghosts(field_velocity_z, 0, 0, 1);

  int gx, gy, gz;

  field_descr->ghosts(field_density, &gx, &gy, &gz);
  unit_assert(gx==3 && gy==3 && gz==3);
  field_descr->ghosts(field_velocity_x, &gx, &gy, &gz);
  unit_assert(gx==1 && gy==0 && gz==0);
  field_descr->ghosts(field_velocity_y, &gx, &gy, &gz);
  unit_assert(gx==0 && gy==1 && gz==0);
  field_descr->ghosts(field_velocity_z, &gx, &gy, &gz);
  unit_assert(gx==0 && gy==0 && gz==1);


  // Minimum value and action


  field_descr->set_minimum (field_density,    1.0, field_action_error);
  field_descr->set_minimum (field_velocity_x, -100.0, field_action_warning);
  field_descr->set_minimum (field_velocity_y, -200.0, field_action_method);
  field_descr->set_minimum (field_velocity_z, -300.0, field_action_timestep);

  field_descr->set_maximum (field_density,    2.0, field_action_error);
  field_descr->set_maximum (field_velocity_x, 100.0, field_action_warning);
  field_descr->set_maximum (field_velocity_y, 200.0, field_action_method);
  field_descr->set_maximum (field_velocity_z, 300.0, field_action_timestep);

  unit_func("min_value");

  unit_assert(field_descr->minimum_value  (field_density)    == 1.0);
  unit_assert(field_descr->minimum_action (field_density)    == field_action_error);
  unit_assert(field_descr->minimum_value  (field_velocity_x) == -100.0);
  unit_assert(field_descr->minimum_action (field_velocity_x) == field_action_warning);
  unit_assert(field_descr->minimum_value  (field_velocity_y) == -200.0);
  unit_assert(field_descr->minimum_action (field_velocity_y) == field_action_method);
  unit_assert(field_descr->minimum_value  (field_velocity_z) == -300.0);
  unit_assert(field_descr->minimum_action (field_velocity_z) == field_action_timestep);

  unit_func("max_value");

  unit_assert(field_descr->maximum_value  (field_density)    == 2.0);
  unit_assert(field_descr->maximum_action (field_density)    == field_action_error);
  unit_assert(field_descr->maximum_value  (field_velocity_x) == 100.0);
  unit_assert(field_descr->maximum_action (field_velocity_x) == field_action_warning);
  unit_assert(field_descr->maximum_value  (field_velocity_y) == 200.0);
  unit_assert(field_descr->maximum_action (field_velocity_y) == field_action_method);
  unit_assert(field_descr->maximum_value  (field_velocity_z) == 300.0);
  unit_assert(field_descr->maximum_action (field_velocity_z) == field_action_timestep);

  unit_finalize();
}
