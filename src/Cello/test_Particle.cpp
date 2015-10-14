// See LICENSE_CELLO file for license and copyright information

/// @file     test_Particle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-10-12
/// @brief    Test program for the Particle class

#include "main.hpp"
#include "test.hpp"

#include "data.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ParticleDescr");
  ParticleDescr particle_descr;

  unit_class("ParticleData");
  ParticleData particle_data;

  Particle particle (&particle_descr, &particle_data);

  //--------------------------------------------------
  //   Type
  //--------------------------------------------------

  // new_type()
  // num_types()

  unit_func ("num_types()");
  unit_assert (particle.num_types() == 0);
  const int i_dark  = particle.new_type ("dark");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 1);
  unit_func ("num_types()");
  unit_assert (i_dark == 0);
  const int i_trace = particle.new_type ("trace");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 2);
  unit_func ("num_types()");
  unit_assert (i_trace == 1);
  const int i_star  = particle.new_type ("star");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 3);
  unit_func ("num_types()");
  unit_assert (i_star == 2);
  const int i_sink  = particle.new_type ("sink");
  unit_func ("new_type()");
  unit_assert (particle.num_types() == 4);
  unit_func ("num_types()");
  unit_assert (i_sink == 3);

  unit_func ("type_name()");
  unit_assert(particle.type_name(i_dark) == "dark");
  unit_assert(particle.type_name(i_trace) == "trace");
  unit_assert(particle.type_name(i_star) == "star");
  unit_assert(particle.type_name(i_sink) == "sink");

  unit_func ("type_index()");
  unit_assert(particle.type_index("dark")  == i_dark);
  unit_assert(particle.type_index("trace") == i_trace);
  unit_assert(particle.type_index("star")  == i_star);
  unit_assert(particle.type_index("sink")  == i_sink);

  //--------------------------------------------------
  //  Attribute
  //--------------------------------------------------

  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 0);
  unit_func ("new_attribute");
  const int id_x = particle.new_attribute (i_dark, "position_x", type_float);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 1);
  unit_func ("new_attribute");
  const int id_y = particle.new_attribute (i_dark, "position_y", type_float);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 2);
  unit_func ("new_attribute");
  const int id_z = particle.new_attribute (i_dark, "position_z", type_float);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 3);
  unit_func ("new_attribute");
  const int id_vx = particle.new_attribute (i_dark, "velocity_x", type_double);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 4);
  unit_func ("new_attribute");
  const int id_vy = particle.new_attribute (i_dark, "velocity_y", type_double);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 5);
  unit_func ("new_attribute");
  const int id_vz = particle.new_attribute (i_dark, "velocity_z", type_double);
  unit_func ("num_attributes");
  unit_assert (particle.num_attributes(i_dark) == 6);
  unit_func ("new_attribute");
  const int id_m = particle.new_attribute (i_dark, "mass",       type_double);
  
  unit_func ("attribute_name()");
  unit_assert(particle.attribute_name(i_dark,id_x) == "position_x");
  unit_assert(particle.attribute_name(i_dark,id_y) == "position_y");
  unit_assert(particle.attribute_name(i_dark,id_z) == "position_z");
  unit_assert(particle.attribute_name(i_dark,id_vx) == "velocity_x");
  unit_assert(particle.attribute_name(i_dark,id_vy) == "velocity_y");
  unit_assert(particle.attribute_name(i_dark,id_vz) == "velocity_z");

  unit_func ("attribute_index()");
  unit_assert(particle.attribute_index(i_dark,"position_x")  == id_x);

  particle.new_attribute (i_trace, "position_x", type_long_double);
  particle.new_attribute (i_trace, "position_y", type_long_double);
  particle.new_attribute (i_trace, "position_z", type_long_double);

  //--------------------------------------------------
  //   Batch
  //--------------------------------------------------

  //--------------------------------------------------
  //   Particle
  //--------------------------------------------------

  //--------------------------------------------------
  //   Grouping
  //--------------------------------------------------


  unit_finalize();

  exit_();


}


PARALLEL_MAIN_END

