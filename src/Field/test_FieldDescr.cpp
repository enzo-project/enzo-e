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
  unit_class ("Field");
  unit_func ("Field");
  FieldDescr density    ("density");
  FieldDescr velocity_x ("velocity_x");
  FieldDescr velocity_y ("velocity_y");
  unit_assert(true);

  unit_func ("name");
  unit_assert(density.name() == "density");
  unit_assert(velocity_x.name() == "velocity_x");
  unit_assert(velocity_y.name() == "velocity_y");

  unit_func ("centering");
  velocity_x.set_centering(1,false);  // Is this the right way around?
  velocity_y.set_centering(0,false);

  unit_assert (density.centering()[0] == true && 
	       density.centering()[1] == true );
  unit_assert (velocity_x.centering()[0] == true && 
	       velocity_x.centering()[1] == false);
  unit_assert (velocity_y.centering()[0] == false && 
	       velocity_y.centering()[1] == true);

  unit_func ("min_value");
  density.set_min_value (1.0);
  unit_assert (density.min_value () == 1.0);

  unit_func ("max_value");
  velocity_x.set_max_value (1.0e6);
  velocity_y.set_max_value (1.0e6);
  unit_assert (velocity_x.max_value () == 1.0e6);
  unit_assert (velocity_y.max_value () == 1.0e6);

  unit_func("min_action");
  density.set_min_action(field_action_assign);
  unit_assert (density.min_action() == field_action_assign);

  unit_func("max_action");
  velocity_x.set_max_action(field_action_warning);
  velocity_y.set_max_action(field_action_error);
  unit_assert (velocity_x.max_action() == field_action_warning);
  unit_assert (velocity_y.max_action() == field_action_error);

  unit_func("precision");
  density.set_precision(precision_default);
  velocity_x.set_precision(precision_single);
  velocity_y.set_precision(precision_double);
  int default_precision = precision_unknown;
#ifdef CONFIG_PRECISION_SINGLE
  default_precision = precision_single;
#endif
#ifdef CONFIG_PRECISION_DOUBLE
  default_precision = precision_double;
#endif
  unit_assert (density.precision()    == default_precision);
  unit_assert (velocity_x.precision() == precision_single);
  unit_assert (velocity_y.precision() == precision_double);

  unit_finalize();
}
