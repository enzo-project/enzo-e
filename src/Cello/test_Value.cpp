// See LICENSE_CELLO file for license and copyright information

/// @file     test_Value.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Value class

#include "main.hpp"
#include "test.hpp"

#include "problem.hpp"

void generate_input()
{
  std::fstream fp;

  fp.open ("test.in",std::fstream::out);

  Group { 
    value1 = [1.0*x+2.0*y-5.0*z];

    value2 = [1.0, (x+y+z > 0.0), 2.0];

    value3 = [1.0, "input/testValue.png"];  // 16 x 8: (x < y) || (x >= y+2.0)
    value4 = [scalar_expr_1, mask_expr_1,
	      scalar_expr_2, mask_png_2];
  }
  fp << "Group {\n";
  fp << "  value_png  = [1.0,\"input/Cello.png\"];\n";
  fp << "  value_x_lt_y = [2.0,x + 0.5*t < 2.0*y - z];\n";
  fp << "}\n";

  fp.close();

}

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Value");

  Value * value = new Value;

  unit_assert (value != NULL);

  //--------------------------------------------------

  unit_func ("function()");

  unit_assert (false);

  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

