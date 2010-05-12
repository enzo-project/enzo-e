// $Id: test_field.cpp 1412 2010-05-05 00:10:09Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_field.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Test program for the FieldDescr class

#include "cello.h"

#include "test.hpp"
#include "field.hpp"

int main()
{

  unit_class ("Field");
  unit_func ("Field");
  FieldDescr * field = new FieldDescr(2);
  unit_assert(field != 0);

  unit_func ("add_field");
  unit_assert(field->field_count() == 0);
  int index_density    = field->add_field("density");
  unit_assert(field->field_count() == 1);
  int index_velocity_x = field->add_field("velocity_x");
  unit_assert(field->field_count() == 2);
  int index_velocity_y = field->add_field("velocity_y");
  unit_assert(field->field_count() == 3);

  unit_func ("name");
  unit_assert(field->name(index_density)    == "density");
  unit_assert(field->name(index_velocity_x) == "velocity_x");
  unit_assert(field->name(index_velocity_y) == "velocity_y");

  unit_func ("index");
  unit_assert(field->index("density") == index_density);
  unit_assert(field->index("velocity_x") == index_velocity_x);
  unit_assert(field->index("velocity_y") == index_velocity_y);

  unit_func ("centering");
  field->set_centering(index_velocity_x,1,false);
  field->set_centering(index_velocity_y,0,false);  // Is this the right way around?

  unit_assert (field->centering(index_density)[0] && 
	       field->centering(index_density)[1]);
  unit_assert (field->centering(index_velocity_x)[0] && 
	       ! field->centering(index_velocity_x)[1]);
  unit_assert (! field->centering(index_velocity_y)[0] && 
	       field->centering(index_velocity_y)[1]);

  unit_func ("min_value");
  field->set_min_value (index_density,1.0);
  unit_assert (field->min_value (index_density) == 1.0);

  unit_func ("max_value");
  field->set_max_value (index_velocity_x,1.0e6);
  field->set_max_value (index_velocity_y,1.0e6);
  unit_assert (field->max_value (index_velocity_x) == 1.0e6);
  unit_assert (field->max_value (index_velocity_y) == 1.0e6);

  unit_func("min_action");
  field->set_min_action(index_density,field_action_set);
  unit_assert (field->min_action(index_density) == field_action_set);

  unit_func("max_action");
  field->set_max_action(index_velocity_x,field_action_warning);
  field->set_max_action(index_velocity_y,field_action_error);
  unit_assert (field->max_action(index_velocity_x) == field_action_warning);
  unit_assert (field->max_action(index_velocity_y) == field_action_error);

  unit_func("precision");
  field->set_precision(index_density,   precision_default);
  field->set_precision(index_velocity_x,precision_32bit);
  field->set_precision(index_velocity_y,precision_64bit);
  int default_precision = precision_unknown;
#ifdef CONFIG_PRECISION_SINGLE
  default_precision = precision_32bit;
#endif
#ifdef CONFIG_PRECISION_DOUBLE
  default_precision = precision_64bit;
#endif
  unit_assert (field->precision(index_density)    == default_precision);
  unit_assert (field->precision(index_velocity_x) == precision_32bit);
  unit_assert (field->precision(index_velocity_y) == precision_64bit);

}
