// See LICENSE_CELLO file for license and copyright information

/// @file     test_Mask.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-29
/// @brief    Test program for the Mask class

#include <fstream>

#include "main.hpp"
#include "test.hpp"

#include "problem.hpp"

void generate_input()
{
  std::fstream fp;

  fp.open ("test.in",std::fstream::out);

  fp << "Group {\n";
  fp << "  value_png  = [1.0,\"input/Cello.png\"];\n";
  fp << "  value_x_lt_y = [2.0,x + 0.5*t < 2.0*y - z];\n";
  fp << "}\n";

  fp.close();

}

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  Parameters parameters;

  generate_input();
  parameters.read("test.in");

  unit_init(0,1);

  unit_class("MaskExpr");

  //--------------------------------------------------

  const int nx = 4;
  const int ny = 8;
  const int nz = 12;
  const double xm = -1.0;
  const double ym = -2.0;
  const double zm = -3.0;
  const double xp =  1.0;
  const double yp =  2.0;
  const double zp =  3.0;
  const double t =   7.0;
  Block block(nx,ny,nz, 1,  xm,xp,ym,yp,zm,zp);

  double x[nx], y[ny], z[nz];

  Mask * mask = new MaskExpr (&parameters, "Group:value_x_lt_y",1);

  unit_func ("Mask()");
  unit_assert (mask != NULL);

  block.field_cells(x,y,z);

  unit_finalize();

  unit_func ("evaluate(ix,iy,iz)");

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	unit_assert (mask->evaluate(t,x[ix],y[iy],z[iz]) == 
		     (x[ix] +0.5*t < 2*y[iy]-z[iz]));
      }
    }
  }

  bool bitmask[nx*ny*nz];

  mask->evaluate(bitmask,t,nx,nx,x,ny,ny,y,nz,nz,z);

  unit_func ("evaluate(mask,nx,ny,nz)");

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix+nx*(iy+ny*iz);
	unit_assert (bitmask[i] ==  (x[ix] +0.5*t < 2*y[iy]-z[iz]));
      }
    }
  }

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

