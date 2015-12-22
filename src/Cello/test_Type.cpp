// See LICENSE_CELLO file for license and copyright information

/// @file     test_Type.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for type information

#include "main.hpp"
#include "test.hpp"

#include "cello.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Type");

  unit_assert (cello::type_bytes[type_single] == sizeof(float));
  unit_assert (cello::type_bytes[type_double] == sizeof(double));
  unit_assert (8*cello::type_bytes[type_extended80] == 80);
  unit_assert (8*cello::type_bytes[type_extended96] == 96);
  unit_assert (8*cello::type_bytes[type_quadruple] == 128);
  unit_assert (cello::type_bytes[type_int8] == sizeof(int8_t));
  unit_assert (cello::type_bytes[type_int16] == sizeof(int16_t));
  unit_assert (cello::type_bytes[type_int32] == sizeof(int32_t));
  unit_assert (cello::type_bytes[type_int64] == sizeof(int64_t));

  unit_assert (strcmp(cello::type_name[type_single],"single") == 0);
  unit_assert (strcmp(cello::type_name[type_double],"double") == 0);
  unit_assert (strcmp(cello::type_name[type_extended80],"extended80") == 0);
  unit_assert (strcmp(cello::type_name[type_extended96],"extended96") == 0);
  unit_assert (strcmp(cello::type_name[type_quadruple],"quadruple") == 0);
  unit_assert (strcmp(cello::type_name[type_int8],"int8") == 0);
  unit_assert (strcmp(cello::type_name[type_int16],"int16") == 0);
  unit_assert (strcmp(cello::type_name[type_int32],"int32") == 0);
  unit_assert (strcmp(cello::type_name[type_int64],"int64") == 0);


  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

