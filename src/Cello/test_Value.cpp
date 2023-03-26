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
double expected_val1(double t, double x, double y, double z) noexcept
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

double expected_val2(double t, double x, double y, double z) noexcept
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

double expected_val3(double t, double x, double y, double z) noexcept
{ return (MASK3_VAL1 ? (EXPR3_VAL1) : (EXPR3_VAL2)); }

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
  fp << "}\n";

  fp.close();

}

//----------------------------------------------------------------------

void test_value_obj_(const Value* value, int num) {

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

  FieldDescr * field_descr = new FieldDescr;
  ParticleDescr * particle_descr = new ParticleDescr;

  Data data(nx,ny,nz, 1,  xm,xp,ym,yp,zm,zp,
	    field_descr, particle_descr);

  double xv[nx], yv[ny], zv[nz];
  double dvalues[n];

  for (int i=0; i<n; i++) dvalues[i] = -999.0;

  data.field_cells(xv,yv,zv);

  unit_assert (value != NULL);

  const char* scalar_fn_name;
  const char* array_fn_name;
  double (*get_expected)(double, double, double, double);

  if (num == 1) {
    scalar_fn_name = "evaluate(scalar) [expr] ";
    array_fn_name  = "evaluate(array) [expr]";
    get_expected = &expected_val1;
  } else if (num == 2) {
    scalar_fn_name = "evaluate(scalar) [expr,mask,expr,mask,expr]";
    array_fn_name  = "evaluate(array) [expr,mask,expr,mask,expr]";
    get_expected = &expected_val2;
  } else if (num == 3) {
    scalar_fn_name = "evaluate(scalar) [expr,mask(png),expr]";
    array_fn_name  = "evaluate(array) [expr,mask(png),expr]";
    get_expected = &expected_val3;
  } else {
    ERROR("test_value_obj_", "the num argument must be 1, 2, or 3");
  }

  // test scalar evaluation
  unit_func (scalar_fn_name);
  unit_assert(value->evaluate(t,xv[0],yv[0],zv[0]) ==
              get_expected(t,xv[0],yv[0],zv[0]));

  // test evaluation over a full array
  unit_func (array_fn_name);
  value->evaluate(dvalues,t,nx,nx,xv,ny,ny,yv,nz,nz,zv);

  bool l_equal = true;
  for (int ix=0; ix<nx; ix++) {
    double x=xv[ix];
    for (int iy=0; iy<ny; iy++) {
      double y=yv[iy];
      for (int iz=0; iz<nz; iz++) {
        int i=ix+nx*(iy+ny*iz);
        double z=zv[iz];
        l_equal = l_equal &&  (dvalues[i] == get_expected(t,x,y,z));
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

    // test case 1
    {
      Value value = Value(&parameters, "Group:value1");
      test_value_obj_(&value, 1);
    }

    //--------------------------------------------------

    // test case 2

    {
      Value value = Value(&parameters, "Group:value2");
      test_value_obj_(&value, 2);
    }

    //--------------------------------------------------

    // test case 3
    {
      Value value = Value(&parameters, "Group:value3");
      test_value_obj_(&value, 3);
    }

    //----------------------------------------------------------------------

  }

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

