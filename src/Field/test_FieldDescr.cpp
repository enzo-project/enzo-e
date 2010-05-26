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
  field_descr->insert_group("centered");
  unit_assert(field_descr->group_count()==1);
  field_descr->insert_group("color");
  unit_func("group_count");
  unit_assert(field_descr->group_count()==2);

  unit_func("group_id");

  int group_centered = field_descr->group_id("centered");
  int group_color    = field_descr->group_id("color");

  unit_assert(field_descr->group_id("centered") == group_centered);
  unit_assert(field_descr->group_id("color")    == group_color);

  unit_func("group_name");

  unit_assert(field_descr->group_name(group_centered) == "centered");
  unit_assert(field_descr->group_name(group_color)    == "color");

  // Fields and groups

  unit_func("set_field_in_group");
  field_descr->set_field_in_group(field_density,group_centered);
  
  

  // void insert_field(std::string name_field) throw();
  // void insert_group(std::string name_group) throw();
  // void set_field_in_group(int id_field, int id_group) throw();
  // void set_alignment(int alignment) throw();
  // void set_padding(int padding) throw();
  // void set_courant(double courant) throw();
  // void set_precision(int id_field, precision_type precision) throw();
  // void set_centering(int id_field, bool cx, bool cy, bool cz) throw();
  // void set_ghosts(int id_field, int gx, int gy, int gz) throw();
  // void set_minimum (int id_field, double min_value, field_action min_action) throw();
  // void set_maximum (int id_field, double max_value, field_action max_action) throw();

  // FieldDescr(const FieldDescr & field_descr) throw();
  // FieldDescr & operator= (const FieldDescr & field_descr) throw();
  // int field_count() const throw();
  // std::string field_name(size_t id_field) const throw(std::out_of_range);
  // int field_id(const std::string name) const throw(std::out_of_range);
  // int group_count() const throw();
  // std::string group_name(int id_group) const throw(std::out_of_range);
  // int group_id(const std::string name) const throw(std::out_of_range);
  // bool field_in_group(int id_field, int id_group) const throw(std::out_of_range);
  // int alignment() const throw();
  // int padding() const throw();
  // void centering(int id_field, bool * cx, bool * cy, bool * cz) const throw(std::out_of_range);
  // void ghosts(int id_field, int * gx, int * gy, int * gz) const throw(std::out_of_range);
  // precision_type precision(int id_field) const throw(std::out_of_range);

  //======================================================================

//   FieldDescr density    ("density");
//   FieldDescr velocity_x ("velocity_x");
//   FieldDescr velocity_y ("velocity_y");
//   unit_assert(true);

//   unit_func ("name");
//   unit_assert(density.name() == "density");
//   unit_assert(velocity_x.name() == "velocity_x");
//   unit_assert(velocity_y.name() == "velocity_y");

//   unit_func ("centering");
//   velocity_x.set_centering(1,false);  // Is this the right way around?
//   velocity_y.set_centering(0,false);

//   unit_assert (density.centering()[0] == true && 
// 	       density.centering()[1] == true );
//   unit_assert (velocity_x.centering()[0] == true && 
// 	       velocity_x.centering()[1] == false);
//   unit_assert (velocity_y.centering()[0] == false && 
// 	       velocity_y.centering()[1] == true);

//   unit_func ("min_value");
//   density.set_min_value (1.0);
//   unit_assert (density.min_value () == 1.0);

//   unit_func ("max_value");
//   velocity_x.set_max_value (1.0e6);
//   velocity_y.set_max_value (1.0e6);
//   unit_assert (velocity_x.max_value () == 1.0e6);
//   unit_assert (velocity_y.max_value () == 1.0e6);

//   unit_func("min_action");
//   density.set_min_action(field_action_assign);
//   unit_assert (density.min_action() == field_action_assign);

//   unit_func("max_action");
//   velocity_x.set_max_action(field_action_warning);
//   velocity_y.set_max_action(field_action_error);
//   unit_assert (velocity_x.max_action() == field_action_warning);
//   unit_assert (velocity_y.max_action() == field_action_error);

//   unit_func("precision");
//   density.set_precision(precision_default);
//   velocity_x.set_precision(precision_single);
//   velocity_y.set_precision(precision_double);
//   int default_precision = precision_unknown;
// #ifdef CONFIG_PRECISION_SINGLE
//   default_precision = precision_single;
// #endif
// #ifdef CONFIG_PRECISION_DOUBLE
//   default_precision = precision_double;
// #endif
//   unit_assert (density.precision()    == default_precision);
//   unit_assert (velocity_x.precision() == precision_single);
//   unit_assert (velocity_y.precision() == precision_double);

  unit_finalize();
}
