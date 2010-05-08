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

  //----------------------------------------------------------------------
  // Particle
  //----------------------------------------------------------------------

  ParticleDescr p0, p1, p2;

  p0.set_name("particle_0");
  p1.set_name("particle_1");
  p2.set_name("particle_2");

  unit_func("add_particle");
  data.add_particle(&p0);
  data.add_particle(&p1);
  data.add_particle(&p2);
  unit_assert(true);

  unit_func("particle_count");
  unit_assert(data.particle_count() == 3);

  unit_func("particle(int)");
  unit_assert(data.particle(0) == &p0);
  unit_assert(data.particle(1) == &p1);
  unit_assert(data.particle(2) == &p2);
  unit_assert(data.particle(3) == 0);
  
  unit_func("particle(string)");
  unit_assert(data.particle("particle_0") == &p0);
  unit_assert(data.particle("particle_1") == &p1);
  unit_assert(data.particle("particle_2") == &p2);
  unit_assert(data.particle("nonexistent particle") == 0);

  //----------------------------------------------------------------------
  // Fields
  //----------------------------------------------------------------------

  FieldDescr f0, f1, f2;

  f0.set_name("field_0");
  f1.set_name("field_1");
  f2.set_name("field_2");

  unit_func("add_field");
  data.add_field(&f0);
  data.add_field(&f1);
  data.add_field(&f2);
  unit_assert(true);

  unit_func("field_count");
  unit_assert(data.field_count() == 3);

  unit_func("field(int)");
  unit_assert(data.field(0) == &f0);
  unit_assert(data.field(1) == &f1);
  unit_assert(data.field(2) == &f2);
  unit_assert(data.field(3) == 0);
  
  unit_func("field(string)");
  unit_assert(data.field("field_0") == &f0);
  unit_assert(data.field("field_1") == &f1);
  unit_assert(data.field("field_2") == &f2);
  unit_assert(data.field("nonexistent field") == 0);
}
