// See LICENSE_CELLO file for license and copyright information

/// @file     test_ScalarData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the ScalarData class

#include "main.hpp"
#include "test.hpp"

#include "data.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ScalarData");

#ifdef NEW_SYNC      

  ScalarData<double> scalar_data;
  ScalarDescr        scalar_descr;
  Scalar<double> scalar (&scalar_descr,&scalar_data);

  //--------------------------------------------------

  unit_func ("new_value()");

  int i0 = scalar.new_value("i0");
  int i1 = scalar.new_value("i1");
  int i2 = scalar.new_value("i2");
  int i3 = scalar.new_value("i3");
  int i4 = scalar.new_value("i4");

  scalar.value(i1) = 1.0;
  scalar.value(i2) = 2.0;
  scalar.value(i3) = 3.0;
  scalar.value(i4) = 4.0;

  unit_assert(scalar.size() == 5);

  unit_func ("value()");
  unit_assert(scalar.value(i0) == 0.0);
  unit_assert(scalar.value(i1) == 1.0);
  unit_assert(scalar.value(i2) == 2.0);
  unit_assert(scalar.value(i3) == 3.0);
  unit_assert(scalar.value(i4) == 4.0);

  unit_func ("name()");
  unit_assert(scalar.name(i0) == "i0");
  unit_assert(scalar.name(i1) == "i1");
  unit_assert(scalar.name(i2) == "i2");
  unit_assert(scalar.name(i3) == "i3");
  unit_assert(scalar.name(i4) == "i4");

  unit_func ("index()");
  unit_assert(scalar.index("i0") == i0);
  unit_assert(scalar.index("i1") == i1);
  unit_assert(scalar.index("i2") == i2);
  unit_assert(scalar.index("i3") == i3);
  unit_assert(scalar.index("i4") == i4);

  unit_func ("assignment");
  scalar.value(i1) = 2.0;
  unit_assert(scalar.value(i1) == 2.0);
  scalar.value(i3) = -3.0;
  unit_assert(scalar.value(i3) == -3.0);
  scalar.value(i3) = scalar.value(i4);
  unit_assert(scalar.value(i3) == scalar.value(i4));
  
  unit_assert(scalar.value(i0) == 0.0);
  unit_assert(scalar.value(i1) == 2.0);
  unit_assert(scalar.value(i2) == 2.0);
  unit_assert(scalar.value(i3) == 4.0);
  unit_assert(scalar.value(i4) == 4.0);

  unit_func ("index() (out of range)");
  unit_assert(scalar.index("i6") == -1);
  //--------------------------------------------------

  unit_finalize();
#endif
  exit_();
  
}

PARALLEL_MAIN_END

