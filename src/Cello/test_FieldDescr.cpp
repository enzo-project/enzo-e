// See LICENSE_CELLO file for license and copyright information

/// @file     test_FieldDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Test program for the FieldDescr class

#include "main.hpp" 
#include "test.hpp"

#include "data.hpp"

struct field_info_type {
  int field_density;
  int field_velocity_x;
  int field_velocity_y;
  int field_velocity_z;
  int field_total_energy;

  int gx, gy, gz;
  int cx, cy, cz;
};

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  //----------------------------------------------------------------------
  unit_init(0,1);
  //----------------------------------------------------------------------

  unit_class("FieldDescr");
  FieldDescr * fd = 0;
  struct field_info_type info;

  unit_func ("FieldDescr");
  fd = new FieldDescr;
  unit_assert(fd != 0);
  printf ("sizeof(FieldDescr) = %lud\n",sizeof(FieldDescr));

  unit_func("~FieldDescr");
  delete fd;
  unit_assert(true);

  fd = new FieldDescr;

  // Fields

  unit_func("insert_permanent");

  int id;

  unit_assert(fd->field_count()==0);

  id = fd->insert_permanent("density");
  unit_func("field_count");
  unit_assert(fd->field_count()==1);
  fd->insert_permanent("velocity_x");
  unit_func("field_count");
  unit_assert(fd->field_count()==2);
  fd->insert_permanent("velocity_y");
  unit_func("field_count");
  unit_assert(fd->field_count()==3);
  fd->insert_permanent("velocity_z");
  unit_func("field_count");
  unit_assert(fd->field_count()==4);

  unit_func("is_permanent");
  id = fd->insert_permanent("total_energy");
  unit_assert(fd->is_permanent(id));
  unit_func("field_count");
  unit_assert(fd->field_count()==5);

  unit_func("is_permanent");
  id = fd->insert_temporary("temporary_1");
  unit_assert(! fd->is_permanent(id));
  unit_func("field_count");
  unit_assert(fd->field_count()==6);
  id = fd->insert_temporary("temporary_2");
  unit_assert(! fd->is_permanent(id));
  unit_func("field_count");
  unit_assert(fd->field_count()==7);

  


  unit_func("field_id");

  info.field_density      = fd->field_id("density");
  info.field_velocity_x   = fd->field_id("velocity_x");
  info.field_velocity_y   = fd->field_id("velocity_y");
  info.field_velocity_z   = fd->field_id("velocity_z");
  info.field_total_energy = fd->field_id("total_energy");

  unit_assert(fd->field_id("density")      == info.field_density);
  unit_assert(fd->field_id("velocity_x")   == info.field_velocity_x);
  unit_assert(fd->field_id("velocity_y")   == info.field_velocity_y);
  unit_assert(fd->field_id("velocity_z")   == info.field_velocity_z);
  unit_assert(fd->field_id("total_energy") == info.field_total_energy);

  unit_func("is_field");

  unit_assert(fd->is_field("density"));
  unit_assert(! fd->is_field("not_a_field"));

  unit_func("field_name");

  unit_assert(fd->field_name(info.field_density)      == "density");
  unit_assert(fd->field_name(info.field_velocity_x)   == "velocity_x");
  unit_assert(fd->field_name(info.field_velocity_y)   == "velocity_y");
  unit_assert(fd->field_name(info.field_velocity_z)   == "velocity_z");
  unit_assert(fd->field_name(info.field_total_energy) == "total_energy");

  //----------------------------------------------------------------------
  // Global attributes
  //----------------------------------------------------------------------

  // (set and reset in case test value is a default)

  fd->set_alignment(8);
  fd->set_padding(64);
  
  unit_func("alignment");
  unit_assert(fd->alignment() == 8);
  unit_func("padding");
  unit_assert(fd->padding() == 64);

  fd->set_alignment(4);
  fd->set_padding(32);
  
  unit_func("alignment");
  unit_assert(fd->alignment() == 4);
  unit_func("padding");
  unit_assert(fd->padding() == 32);
  
  //----------------------------------------------------------------------
  // Field attributes
  //----------------------------------------------------------------------

  // Precision

  unit_func("precision");

  fd->set_precision(info.field_density,    precision_single);
  fd->set_precision(info.field_velocity_x, precision_double);
  fd->set_precision(info.field_velocity_y, precision_double);
  fd->set_precision(info.field_velocity_z, precision_double);

  unit_assert(fd->precision(info.field_density)      == precision_single);
  unit_assert(fd->precision(info.field_velocity_x)   == precision_double);
  unit_assert(fd->precision(info.field_velocity_y)   == precision_double);
  unit_assert(fd->precision(info.field_velocity_z)   == precision_double);
  unit_assert(fd->precision(info.field_total_energy) == default_precision);
  

  unit_func("bytes_per_element");
  unit_assert(fd->bytes_per_element(info.field_density)==4);

  // Centering

  unit_func("centering");

  fd->set_centering(info.field_velocity_x, 1, 0, 0);
  fd->set_centering(info.field_velocity_y, 0, 1, 0);
  fd->set_centering(info.field_velocity_z, 0, 0, 1);


  fd->centering(info.field_density, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==0 && info.cy==0 && info.cz==0);

  fd->centering(info.field_velocity_x, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==1 && info.cy==0 && info.cz==0);

  fd->centering(info.field_velocity_y, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==0 && info.cy==1 && info.cz==0);

  fd->centering(info.field_velocity_z, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==0 && info.cy==0 && info.cz==1);
  
  // Ghost zone depth

  unit_func("ghosts");

  fd->set_ghost_depth(info.field_density, 3, 3, 3);
  fd->set_ghost_depth(info.field_velocity_x, 1, 0, 0);
  fd->set_ghost_depth(info.field_velocity_y, 0, 1, 0);
  fd->set_ghost_depth(info.field_velocity_z, 0, 0, 1);

  fd->ghost_depth(info.field_density, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==3 && info.gy==3 && info.gz==3);
  fd->ghost_depth(info.field_velocity_x, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==1 && info.gy==0 && info.gz==0);
  fd->ghost_depth(info.field_velocity_y, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==1 && info.gz==0);
  fd->ghost_depth(info.field_velocity_z, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==0 && info.gz==1);


  //======================================================================
  // BIG THREE
  //======================================================================

  // Assign
  FieldDescr fd_assign;
  fd_assign = *fd;

  // Copy
  FieldDescr fd_copy (*fd);

  // Delete original to check for deep copy

  unit_func("~FieldDescr");
  delete fd;
  fd = 0;

  unit_func("assign:field_count");
  unit_assert(fd_assign.field_count()==7);

  unit_func("assign:field_id");
  unit_assert(fd_assign.field_id("density")      == info.field_density);
  unit_assert(fd_assign.field_id("velocity_x")   == info.field_velocity_x);
  unit_assert(fd_assign.field_id("velocity_y")   == info.field_velocity_y);
  unit_assert(fd_assign.field_id("velocity_z")   == info.field_velocity_z);
  unit_assert(fd_assign.field_id("total_energy") == info.field_total_energy);

  unit_func("assign:field_name");
  unit_assert(fd_assign.field_name(info.field_density)      == "density");
  unit_assert(fd_assign.field_name(info.field_velocity_x)   == "velocity_x");
  unit_assert(fd_assign.field_name(info.field_velocity_y)   == "velocity_y");
  unit_assert(fd_assign.field_name(info.field_velocity_z)   == "velocity_z");
  unit_assert(fd_assign.field_name(info.field_total_energy) == "total_energy");


  unit_func("assign:alignment");
  unit_assert(fd_assign.alignment() == 4);
  unit_func("assign:padding");
  unit_assert(fd_assign.padding() == 32);

  unit_func("assign:precision");

  unit_assert(fd_assign.precision(info.field_density)      == precision_single);
  unit_assert(fd_assign.precision(info.field_velocity_x)   == precision_double);
  unit_assert(fd_assign.precision(info.field_velocity_y)   == precision_double);
  unit_assert(fd_assign.precision(info.field_velocity_z)   == precision_double);
  unit_assert(fd_assign.precision(info.field_total_energy) == default_precision);

  unit_func("assign:centering");
  fd_assign.centering(info.field_density, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==0 && info.cy==0 && info.cz==0);

  fd_assign.centering(info.field_velocity_x, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==1 && info.cy==0 && info.cz==0);

  fd_assign.centering(info.field_velocity_y, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==0 && info.cy==1 && info.cz==0);

  fd_assign.centering(info.field_velocity_z, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==0 && info.cy==0 && info.cz==1);

  unit_func("assign:ghosts");

  fd_assign.ghost_depth(info.field_density, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==3 && info.gy==3 && info.gz==3);
  fd_assign.ghost_depth(info.field_velocity_x, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==1 && info.gy==0 && info.gz==0);
  fd_assign.ghost_depth(info.field_velocity_y, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==1 && info.gz==0);
  fd_assign.ghost_depth(info.field_velocity_z, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==0 && info.gz==1);


  unit_func("copy:FieldDescr(FieldDescr)");
  unit_assert(fd_copy.field_count()==7);
  printf ("%s:%d num_permanent = %d",__FILE__,__LINE__,fd_copy.num_permanent());
  unit_assert(fd_copy.num_permanent()==5);

  unit_func("copy:field_id");
  unit_assert(fd_copy.field_id("density")      == info.field_density);
  unit_assert(fd_copy.field_id("velocity_x")   == info.field_velocity_x);
  unit_assert(fd_copy.field_id("velocity_y")   == info.field_velocity_y);
  unit_assert(fd_copy.field_id("velocity_z")   == info.field_velocity_z);
  unit_assert(fd_copy.field_id("total_energy") == info.field_total_energy);

  unit_func("copy:field_name");
  unit_assert(fd_copy.field_name(info.field_density)      == "density");
  unit_assert(fd_copy.field_name(info.field_velocity_x)   == "velocity_x");
  unit_assert(fd_copy.field_name(info.field_velocity_y)   == "velocity_y");
  unit_assert(fd_copy.field_name(info.field_velocity_z)   == "velocity_z");
  unit_assert(fd_copy.field_name(info.field_total_energy) == "total_energy");


  unit_func("copy:alignment");
  unit_assert(fd_copy.alignment() == 4);
  unit_func("copy:padding");
  unit_assert(fd_copy.padding() == 32);

  unit_func("copy:precision");

  unit_assert(fd_copy.precision(info.field_density)      == precision_single);
  unit_assert(fd_copy.precision(info.field_velocity_x)   == precision_double);
  unit_assert(fd_copy.precision(info.field_velocity_y)   == precision_double);
  unit_assert(fd_copy.precision(info.field_velocity_z)   == precision_double);
  unit_assert(fd_copy.precision(info.field_total_energy) == default_precision);

  unit_func("copy:centering");

  fd_copy.centering(info.field_density, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==0 && info.cy==0 && info.cz==0);

  fd_copy.centering(info.field_velocity_x, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==1 && info.cy==0 && info.cz==0);

  fd_copy.centering(info.field_velocity_y, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==0 && info.cy==1 && info.cz==0);

  fd_copy.centering(info.field_velocity_z, &info.cx, &info.cy, &info.cz);
  unit_assert(info.cx==0 && info.cy==0 && info.cz==1);

  unit_func("copy:ghosts");

  fd_copy.ghost_depth(info.field_density, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==3 && info.gy==3 && info.gz==3);
  fd_copy.ghost_depth(info.field_velocity_x, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==1 && info.gy==0 && info.gz==0);
  fd_copy.ghost_depth(info.field_velocity_y, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==1 && info.gz==0);
  fd_copy.ghost_depth(info.field_velocity_z, &info.gx, &info.gy, &info.gz);
  unit_assert(info.gx==0 && info.gy==0 && info.gz==1);


  //----------------------------------------------------------------------
  unit_finalize();
  //----------------------------------------------------------------------

  exit_();
}
PARALLEL_MAIN_END
