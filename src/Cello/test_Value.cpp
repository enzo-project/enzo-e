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

#define EXPR3_VAL1  (t + 10.0*x + 100.0*y + 1000.0*z)
#define EXPR3_STR1 "(t + 10.0*x + 100.0*y + 1000.0*z)"
#define MASK3_VAL1  (x + y >= 1.99 || y - x > 2.001)
#define MASK3_STR1  "\"input/testValue.png\""
#define EXPR3_VAL2  (1.0 - t - 10.0*x - 100.0*y - 1000.0*z)
#define EXPR3_STR2 "(1.0 - t - 10.0*x - 100.0*y - 1000.0*z)"
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

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  Parameters parameters;

  generate_input();
  parameters.read("test.in");

  unit_class("Value");


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
  Block block(field_descr,nx,ny,nz, 1,  xm,xp,ym,yp,zm,zp);

  double xv[nx], yv[ny], zv[nz];
  double dvalues[n];

  for (int i=0; i<n; i++) dvalues[i] = -999.0;

  block.field_cells(xv,yv,zv);

  double x=xv[0];
  double y=yv[0];
  double z=zv[0];
    
  unit_func ("evaluate(scalar) [expr] ");

  Value * value1 = new Value(&parameters, "Group:value1");
  unit_assert (value1 != NULL);

  unit_assert(value1->evaluate(t,xv[0],yv[0],zv[0]) == EXPR1_VAL);

  value1->evaluate(dvalues,t,nx,nx,xv,ny,ny,yv,nz,nz,zv);

  unit_func ("evaluate(array) [expr]");

  for (int ix=0; ix<nx; ix++) {
    double x=xv[ix];
    for (int iy=0; iy<ny; iy++) {
      double y=yv[iy];
      for (int iz=0; iz<nz; iz++) {
	int i=ix+nx*(iy+ny*iz);
	double z=zv[iz];
	  unit_assert (dvalues[i] == EXPR1_VAL);
	}
    }
  }
  
  
  //--------------------------------------------------

  unit_func ("evaluate(scalar) [expr,mask,expr] ");

  Value * value2 = new Value(&parameters, "Group:value2");
  unit_assert (value2 != NULL);

  double val2 = (MASK2_VAL1 ? (EXPR2_VAL1) : ( (MASK2_VAL2 ? (EXPR2_VAL2) : (EXPR2_VAL3))));
  unit_assert(value2->evaluate(t,xv[0],yv[0],zv[0]) == val2);

  for (int i=0; i<n; i++) dvalues[i] = -999.0;
  value2->evaluate(dvalues,t,nx,nx,xv,ny,ny,yv,nz,nz,zv);

  unit_func ("evaluate(array) [expr,mask,expr]");

  for (int ix=0; ix<nx; ix++) {
    double x=xv[ix];
    for (int iy=0; iy<ny; iy++) {
      double y=yv[iy];
      for (int iz=0; iz<nz; iz++) {
	int i=ix+nx*(iy+ny*iz);
	double z=zv[iz];
	val2 = (MASK2_VAL1 ? (EXPR2_VAL1) : ( (MASK2_VAL2 ? (EXPR2_VAL2) : (EXPR2_VAL3))));
	unit_assert (dvalues[i] == val2);
      }
    }
  }

  //--------------------------------------------------

  unit_func ("evaluate(scalar) [expr,mask(png),expr] ");

  Value * value3 = new Value(&parameters, "Group:value3");
  unit_assert (value3 != NULL);

  double val3 = (MASK3_VAL1 ? (EXPR3_VAL1) : (EXPR3_VAL2));
  unit_assert(value3->evaluate(t,xv[0],yv[0],zv[0]) == val3);

  for (int i=0; i<n; i++) dvalues[i] = -999.0;
  value3->evaluate(dvalues,t,nx,nx,xv,ny,ny,yv,nz,nz,zv);

  unit_func ("evaluate(array) [expr,mask(png),expr]");

  for (int ix=0; ix<nx; ix++) {
    double x=xv[ix];
    for (int iy=0; iy<ny; iy++) {
      double y=yv[iy];
      for (int iz=0; iz<nz; iz++) {
	int i=ix+nx*(iy+ny*iz);
	double z=zv[iz];
	val3 = (MASK3_VAL1 ? (EXPR3_VAL1) : (EXPR3_VAL2));
	unit_assert (dvalues[i] == val3);
      }
    }
  }

  //----------------------------------------------------------------------

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

