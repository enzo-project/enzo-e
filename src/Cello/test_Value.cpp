// See LICENSE_CELLO file for license and copyright information

/// @file     test_Value.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Value class

#include <fstream>

#include "main.hpp"
#include "test.hpp"

#include "problem.hpp"

#define EXPR1_VAL  (1.0*x + 2.0*y - 5.0*z + t)
#define EXPR1_STR "(1.0*x + 2.0*y - 5.0*z + t)"

double expected1(double t, double x, double y, double z) noexcept
{ return EXPR1_VAL; }

#define EXPR2_VAL1  (1.0 + 2.0*x + 4.0*y + 8.0*z + 16.0*t)
#define EXPR2_STR1 "(1.0 + 2.0*x + 4.0*y + 8.0*z + 16.0*t)"
#define MASK2_VAL1  (x < y)
#define MASK2_STR1 "(x < y)"
#define EXPR2_VAL2  (-10.0 + 1.0*x + 2.0*y + 4.0*z + 8.0*t)
#define EXPR2_STR2 "(-10.0 + 1.0*x + 2.0*y + 4.0*z + 8.0*t)"
#define MASK2_VAL2  (x > 0.0)
#define MASK2_STR2 "(x > 0.0)"
#define EXPR2_VAL3  (100.0 - 1.0*x - 2.0*y - 5.0*z - t)
#define EXPR2_STR3 "(100.0 - 1.0*x - 2.0*y - 5.0*z - t)"

double expected2(double t, double x, double y, double z) noexcept
{
  return (MASK2_VAL1 ? (EXPR2_VAL1)
                     : ( (MASK2_VAL2 ? (EXPR2_VAL2)
                                     : (EXPR2_VAL3) ))
         );
}

#define EXPR3_VAL1  (t + 10.0*x + 100.0*y + 1000.0*z)
#define EXPR3_STR1 "(t + 10.0*x + 100.0*y + 1000.0*z)"
#define MASK3_VAL1  (x + y >= 1.99 || y - x > 2.001)
#define MASK3_STR1  "\"input/testValue.png\""
#define EXPR3_VAL2  (1.0 - t - 10.0*x - 100.0*y - 1000.0*z)
#define EXPR3_STR2 "(1.0 - t - 10.0*x - 100.0*y - 1000.0*z)"

double expected3(double t, double x, double y, double z) noexcept
{ return (MASK3_VAL1 ? (EXPR3_VAL1) : (EXPR3_VAL2)); }

#define EXPR4_VAL (17.0)
#define EXPR4_STR "17.0"

double expected4(double t, double x, double y, double z) noexcept
{ return EXPR4_VAL; }

#define EXPR5_VAL1 (17.0)
#define EXPR5_STR1 "17.0"
#define MASK5_VAL1  (x < y)
#define MASK5_STR1 "(x < y)"
#define EXPR5_VAL2  (1.0 + 2.0*x + 4.0*y + 8.0*z + 16.0*t)
#define EXPR5_STR2 "(1.0 + 2.0*x + 4.0*y + 8.0*z + 16.0*t)"

double expected5(double t, double x, double y, double z) noexcept
{ return MASK5_VAL1 ? (EXPR5_VAL1) : (EXPR5_VAL2); }

//----------------------------------------------------------------------

void generate_input()
{
  std::fstream fp;

  fp.open ("test.in",std::fstream::out);

  fp << "  Domain {\n";
  fp << "     lower = [-4.0, -2.0, -3.0];\n";
  fp << "     upper = [ 4.0,  2.0,  3.0];\n";
  fp << "}\n";

  fp << "  Group {\n";
  fp << "    value1 = [" EXPR1_STR "];  \n";
  fp << "    value2 = [" EXPR2_STR1 ",\n" MASK2_STR1 ",\n" EXPR2_STR2 ",\n" MASK2_STR2 ",\n" EXPR2_STR3 "];\n";
  fp << "    value3 = [" EXPR3_STR1 ",\n" MASK3_STR1 ",\n" EXPR3_STR2 "];\n";
  fp << "    value4 = [" EXPR4_STR "];  \n";
  fp << "    value5 = [" EXPR5_STR1 ",\n" MASK5_STR1 ",\n" EXPR5_STR2 "];\n";
  fp << "}\n";

  fp.close();

}

//----------------------------------------------------------------------

void test_value_obj_(const Value& value, int num, std::string prefix = "") {

  const int nx = 16;
  const int ny = 8;
  const int nz = 12;
  const int n = nx*ny*nz;
  const double xm = -4.0;
  const double ym = -2.0;
  const double zm = -3.0;
  const double xp =  4.0;
  const double yp =  2.0;
  const double zp =  3.0;
  const double t =   7.0;

  // setup xv, yv, zv
  double xv[nx], yv[ny], zv[nz];
  {
    FieldDescr field_descr;
    ParticleDescr particle_descr;

    Data data(nx,ny,nz, 1,  xm,xp,ym,yp,zm,zp,
              &field_descr, &particle_descr);

    data.field_cells(xv,yv,zv);
  }

  // setup array where outputs will get stored
  double dvalues[n];
  for (int i=0; i<n; i++) dvalues[i] = -999.0;

  // get the description and func ptr that yields expected value for test case
  std::string descr;
  double (*expect_fn)(double, double, double, double);
  switch (num) {
    case 1: descr = "[expr]";                     expect_fn = &expected1; break;
    case 2: descr = "[expr,mask,expr,mask,expr]"; expect_fn = &expected2; break;
    case 3: descr = "[expr,mask(png),expr]";      expect_fn = &expected3; break;
    case 4: descr = "[float]";                    expect_fn = &expected4; break;
    case 5: descr = "[float,mask,expr]";          expect_fn = &expected5; break;
    default:  ERROR("test_value_obj_", "the num argument must be 1-5");
  }

  // this function explicitly makes a copy of name to avoid lifetime issues
  auto set_unit_func = [](std::string name) { unit_func(name.c_str()); };

  // now, actually perform the tests:

  // test scalar evaluation
  set_unit_func ( prefix + "evaluate(scalar) " + descr);
  unit_assert(value.evaluate(t,xv[0],yv[0],zv[0]) ==
              expect_fn(t,xv[0],yv[0],zv[0]));

  // test evaluation over a full array
  set_unit_func ( prefix + "evaluate(array) " + descr);
  value.evaluate(dvalues,t,nx,nx,xv,ny,ny,yv,nz,nz,zv);

  bool l_equal = true;
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
        int i=ix+nx*(iy+ny*iz);
        l_equal = l_equal && (dvalues[i] == expect_fn(t,xv[ix],yv[iy],zv[iz]));
      }
    }
  }
  unit_assert (l_equal);

}

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  // we explicitly put functions in a separate scope to ensure that the
  // constructor runs for the Parameters object and each of the Value objects.
  // - This explicitly happens to explicitly try to catch a hypothetical bug at
  //   one time where the parameters object and underlying members in the Value
  //   object may have deleted a Param object more than once (in practice this
  //   never happened - because Value's destructor didn't deallocate the
  //   contained objects)
  {
    Parameters parameters;

    generate_input();
    parameters.read("test.in");

    unit_class("Value");

    // iterate over the test cases:
    for (int i = 1; i <= 5; i++){
      std::string param_str = "Group:value" + std::to_string(i);

      // the curly braces are used to tell the compiler it can destroy the
      // objects outsude of the region
      {
        // first, check the Value object under normal conditions
        Value orig_value = Value(&parameters, param_str);
        test_value_obj_(orig_value, i);

        // next, check a move-constructed Value object
        Value move_constructed(std::move(orig_value));
        test_value_obj_(move_constructed, i, "move-constructed: ");

        // now, let's see what with a move-assigned Value object
        Value move_assigned; // default constructor
        move_assigned = std::move(move_constructed);
        test_value_obj_(move_assigned, i, "move-assigned: ");
      }

      // now, let's explicitly check again with a freshly constructed value
      // object (this should help detect any issues with the destructor)
      {
        // first, check the Value object under normal conditions
        Value value = Value(&parameters, param_str);
        test_value_obj_(value, i, "freshly-reconstructed: ");
      }
    }

  }

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

