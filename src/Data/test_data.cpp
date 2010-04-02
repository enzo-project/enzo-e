#include "cello.h"

#include "error.hpp"
#include "test.hpp"
#include "data.hpp"
#include "particles.hpp"
#include "field.hpp"

int main()
{

  unit_class ("Data");

  unit_open();


  Data data;

  unit_func("set_particles");
  Particles particles;
  data.set_particles(&particles);
  unit_assert(true);

  unit_func("get_particles");
  Particles * p;
  data.get_particles(&p);
  unit_assert(p == & particles);

  unit_func("set_field");
  Field field;
  data.set_fields(&field);
  unit_assert(true);

  unit_func("get_field");
  Field * f;
  data.get_fields(&f);
  unit_assert(f == & field);

  unit_close();

}
