// See LICENSE_CELLO file for license and copyright information

/// @file     test_Type.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for type information

#include "main.hpp"
#include "test.hpp"

#include "cello.hpp"

#include <cstdint>

void check_get_type_enum(){
  // first check floats:
  unit_assert (cello::get_type_enum<float>() == type_single);
  unit_assert (cello::get_type_enum<double>() == type_double);
  {
    int long_double_val = cello::get_type_enum<long double>();
    unit_assert((long_double_val == type_extended80) ||
                (long_double_val == type_extended96) ||
                (long_double_val == type_quadruple));
  }

  // next, check other types:
  unit_assert (cello::get_type_enum<char>() == type_int8);
  unit_assert (cello::get_type_enum<short>() == type_int16);
  unit_assert (cello::get_type_enum<int>() == type_int32);
  unit_assert (cello::get_type_enum<long long>() == type_int64);

  // fixed width integers are considered distinct types from
  // char, short, int, long long (even when they are the same size)
  //
  // so, the following doesn't work at the time of writing this test
  //unit_assert (cello::get_type_enum<std::int8_t>() == type_int8);
  //unit_assert (cello::get_type_enum<std::int16_t>() == type_int16);
  //unit_assert (cello::get_type_enum<std::int32_t>() == type_int32); 
  //unit_assert (cello::get_type_enum<std::int64_t>() == type_int64);

  // finally, check that this works with user defined type aliases
  {
    typedef float my_custom_float1;
    unit_assert (cello::get_type_enum<my_custom_float1>() == type_single);
  }
  {
    using my_custom_float2 = double;
    unit_assert (cello::get_type_enum<my_custom_float2>() == type_double);
  }
}

void check_enum_conversion(){
  // this is mostly just a sanity check
  unit_assert ( type_unknown ==
                cello::convert_enum_precision_to_type(precision_unknown) );
  unit_assert ( type_default ==
                cello::convert_enum_precision_to_type(precision_default) );
  unit_assert ( type_single ==
                cello::convert_enum_precision_to_type(precision_single) );
  unit_assert ( type_double ==
                cello::convert_enum_precision_to_type(precision_double) );
  unit_assert ( type_extended80 ==
                cello::convert_enum_precision_to_type(precision_extended80) );
  unit_assert ( type_extended96 ==
                cello::convert_enum_precision_to_type(precision_extended96) );
  unit_assert ( type_quadruple ==
                cello::convert_enum_precision_to_type(precision_quadruple) );
}

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

  check_get_type_enum();
  check_enum_conversion();

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

