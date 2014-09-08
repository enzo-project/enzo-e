// See LICENSE_CELLO file for license and copyright information

/// @file     test_FieldDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Test program for the FieldDescr class

#include "main.hpp" 
#include "test.hpp"

#include "field.hpp"

struct field_info_type {
  int field_density;
  int field_velocity_x;
  int field_velocity_y;
  int field_velocity_z;
  int field_total_energy;

  int group_density;
  int group_vector;

  int gx, gy, gz;
  bool cx, cy, cz;
};

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  //----------------------------------------------------------------------
  unit_init(0,1);
  //----------------------------------------------------------------------

  unit_class("FieldDescr");
  FieldDescr * field_descr = 0;
  struct field_info_type info;

  unit_func ("FieldDescr");
  field_descr = new FieldDescr;
  unit_assert(field_descr != 0);
  printf ("sizeof(FieldDescr) = %lud\n",sizeof(FieldDescr));

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
  unit_assert(field_descr->field_count()==5);
  field_descr->insert_field("total_energy");
  unit_assert(field_descr->field_count()==5);

  unit_func("field_count");
  unit_assert(field_descr->field_count()==5);

  unit_func("field_id");

  info.field_density      = field_descr->field_id("density");
  info.field_velocity_x   = field_descr->field_id("velocity_x");
  info.field_velocity_y   = field_descr->field_id("velocity_y");
  info.field_velocity_z   = field_descr->field_id("velocity_z");
  info.field_total_energy = field_descr->field_id("total_energy");

  unit_assert(field_descr->field_id("density")      == info.field_density);
  unit_assert(field_descr->field_id("velocity_x")   == info.field_velocity_x);
  unit_assert(field_descr->field_id("velocity_y")   == info.field_velocity_y);
  unit_assert(field_descr->field_id("velocity_z")   == info.field_velocity_z);
  unit_assert(field_descr->field_id("total_energy") == info.field_total_energy);

  unit_func("is_field");

  unit_assert(field_descr->is_field("density"));
  unit_assert(! field_descr->is_field("not_a_field"));

  unit_func("field_name");

  unit_assert(field_descr->field_name(info.field_density)      == "density");
  unit_assert(field_descr->field_name(info.field_velocity_x)   == "velocity_x");
  unit_assert(field_descr->field_name(info.field_velocity_y)   == "velocity_y");
  unit_assert(field_descr->field_name(info.field_velocity_z)   == "velocity_z");
  unit_assert(field_descr->field_name(info.field_total_energy) == "total_energy");

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

  field_descr->set_precision(info.field_density,    precision_single);
  field_descr->set_precision(info.field_velocity_x, precision_double);
  field_descr->set_precision(info.field_velocity_y, precision_double);
  field_descr->set_precision(info.field_velocity_z, precision_double);

  unit_assert(field_descr->precision(info.field_density)      == precision_single);
  unit_assert(field_descr->precision(info.field_velocity_x)   == precision_double);
  unit_assert(field_descr->precision(info.field_velocity_y)   == precision_double);
  unit_assert(field_descr->precision(info.field_velocity_z)   == precision_double);
  unit_assert(field_descr->precision(info.field_total_energy) == default_precision);
  

  unit_func("bytes_per_element");
  unit_assert(field_descr->bytes_per_element(info.field_density)==4);

  // Centering

  unit_func("centering");

  field_descr->set_centering(info.field_velocity_x, false, true, true);
  field_descr->set_centering(info.field_velocity_y, true, false, true);
  field_descr->set_centering(info.field_velocity_z, true,  true, false);


  field_descr->centering(info.field_density, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx && info.cy && info.cz);

  field_descr->centering(info.field_velocity_x, &info.cx, &info.cy, &info.cz);
  unit_assert(! info.cx &&   info.cy &&   info.cz);

  field_descr->centering(info.field_velocity_y, &info.cx, &info.cy, &info.cz);
  unit_assert(  info.cx && ! info.cy &&   info.cz);

  field_descr->centering(info.field_velocity_z, &info.cx, &info.cy, &info.cz);
  unit_assert(  info.cx &&   info.cy && ! info.cz);
  
  // Ghost zone depth

  unit_func("ghosts");

  field_descr->set_ghosts(info.field_density, 3, 3, 3);
  field_descr->set_ghosts(info.field_velocity_x, 1, 0, 0);
  field_descr->set_ghosts(info.field_velocity_y, 0, 1, 0);
  field_descr->set_ghosts(info.field_velocity_z, 0, 0, 1);

  field_descr->ghosts(info.field_density, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==3 && info.gy==3 && info.gz==3);
  field_descr->ghosts(info.field_velocity_x, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==1 && info.gy==0 && info.gz==0);
  field_descr->ghosts(info.field_velocity_y, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==1 && info.gz==0);
  field_descr->ghosts(info.field_velocity_z, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==0 && info.gz==1);


  //======================================================================
  // BIG THREE
  //======================================================================

  // Assign
  FieldDescr field_descr_assign;
  field_descr_assign = *field_descr;

  // Copy
  FieldDescr field_descr_copy (*field_descr);

  // Delete original to check for deep copy

  unit_func("~FieldDescr");
  delete field_descr;
  field_descr = 0;

  unit_func("assign:field_count");
  unit_assert(field_descr_assign.field_count()==5);

  unit_func("assign:field_id");
  unit_assert(field_descr_assign.field_id("density")      == info.field_density);
  unit_assert(field_descr_assign.field_id("velocity_x")   == info.field_velocity_x);
  unit_assert(field_descr_assign.field_id("velocity_y")   == info.field_velocity_y);
  unit_assert(field_descr_assign.field_id("velocity_z")   == info.field_velocity_z);
  unit_assert(field_descr_assign.field_id("total_energy") == info.field_total_energy);

  unit_func("assign:field_name");
  unit_assert(field_descr_assign.field_name(info.field_density)      == "density");
  unit_assert(field_descr_assign.field_name(info.field_velocity_x)   == "velocity_x");
  unit_assert(field_descr_assign.field_name(info.field_velocity_y)   == "velocity_y");
  unit_assert(field_descr_assign.field_name(info.field_velocity_z)   == "velocity_z");
  unit_assert(field_descr_assign.field_name(info.field_total_energy) == "total_energy");


  unit_func("assign:alignment");
  unit_assert(field_descr_assign.alignment() == 4);
  unit_func("assign:padding");
  unit_assert(field_descr_assign.padding() == 32);
  unit_func("assign:courant");
  unit_assert(field_descr_assign.courant() == 0.75);

  unit_func("assign:precision");

  unit_assert(field_descr_assign.precision(info.field_density)      == precision_single);
  unit_assert(field_descr_assign.precision(info.field_velocity_x)   == precision_double);
  unit_assert(field_descr_assign.precision(info.field_velocity_y)   == precision_double);
  unit_assert(field_descr_assign.precision(info.field_velocity_z)   == precision_double);
  unit_assert(field_descr_assign.precision(info.field_total_energy) == default_precision);

  unit_func("assign:centering");
  field_descr_assign.centering(info.field_density, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx && info.cy && info.cz);

  field_descr_assign.centering(info.field_velocity_x, &info.cx, &info.cy, &info.cz);
  unit_assert(! info.cx &&   info.cy &&   info.cz);

  field_descr_assign.centering(info.field_velocity_y, &info.cx, &info.cy, &info.cz);
  unit_assert(  info.cx && ! info.cy &&   info.cz);

  field_descr_assign.centering(info.field_velocity_z, &info.cx, &info.cy, &info.cz);
  unit_assert(  info.cx &&   info.cy && ! info.cz);

  unit_func("assign:ghosts");

  field_descr_assign.ghosts(info.field_density, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==3 && info.gy==3 && info.gz==3);
  field_descr_assign.ghosts(info.field_velocity_x, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==1 && info.gy==0 && info.gz==0);
  field_descr_assign.ghosts(info.field_velocity_y, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==1 && info.gz==0);
  field_descr_assign.ghosts(info.field_velocity_z, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==0 && info.gz==1);


  unit_func("copy:FieldDescr(FieldDescr)");
  unit_assert(field_descr_copy.field_count()==5);

  unit_func("copy:field_id");
  unit_assert(field_descr_copy.field_id("density")      == info.field_density);
  unit_assert(field_descr_copy.field_id("velocity_x")   == info.field_velocity_x);
  unit_assert(field_descr_copy.field_id("velocity_y")   == info.field_velocity_y);
  unit_assert(field_descr_copy.field_id("velocity_z")   == info.field_velocity_z);
  unit_assert(field_descr_copy.field_id("total_energy") == info.field_total_energy);

  unit_func("copy:field_name");
  unit_assert(field_descr_copy.field_name(info.field_density)      == "density");
  unit_assert(field_descr_copy.field_name(info.field_velocity_x)   == "velocity_x");
  unit_assert(field_descr_copy.field_name(info.field_velocity_y)   == "velocity_y");
  unit_assert(field_descr_copy.field_name(info.field_velocity_z)   == "velocity_z");
  unit_assert(field_descr_copy.field_name(info.field_total_energy) == "total_energy");


  unit_func("copy:alignment");
  unit_assert(field_descr_copy.alignment() == 4);
  unit_func("copy:padding");
  unit_assert(field_descr_copy.padding() == 32);
  unit_func("copy:courant");
  unit_assert(field_descr_copy.courant() == 0.75);

  unit_func("copy:precision");

  unit_assert(field_descr_copy.precision(info.field_density)      == precision_single);
  unit_assert(field_descr_copy.precision(info.field_velocity_x)   == precision_double);
  unit_assert(field_descr_copy.precision(info.field_velocity_y)   == precision_double);
  unit_assert(field_descr_copy.precision(info.field_velocity_z)   == precision_double);
  unit_assert(field_descr_copy.precision(info.field_total_energy) == default_precision);

  unit_func("copy:centering");

  field_descr_copy.centering(info.field_density, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx && info.cy && info.cz);

  field_descr_copy.centering(info.field_velocity_x, &info.cx, &info.cy, &info.cz);
  unit_assert(! info.cx &&   info.cy &&   info.cz);

  field_descr_copy.centering(info.field_velocity_y, &info.cx, &info.cy, &info.cz);
  unit_assert(  info.cx && ! info.cy &&   info.cz);

  field_descr_copy.centering(info.field_velocity_z, &info.cx, &info.cy, &info.cz);
  unit_assert(  info.cx &&   info.cy && ! info.cz);

  unit_func("copy:ghosts");

  field_descr_copy.ghosts(info.field_density, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==3 && info.gy==3 && info.gz==3);
  field_descr_copy.ghosts(info.field_velocity_x, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==1 && info.gy==0 && info.gz==0);
  field_descr_copy.ghosts(info.field_velocity_y, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==1 && info.gz==0);
  field_descr_copy.ghosts(info.field_velocity_z, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==0 && info.gz==1);


  //----------------------------------------------------------------------
  unit_finalize();
  //----------------------------------------------------------------------

  exit_();
}
PARALLEL_MAIN_END
