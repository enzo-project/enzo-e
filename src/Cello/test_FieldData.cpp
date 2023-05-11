// See LICENSE_CELLO file for license and copyright information

/// @file     test_FieldData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the FieldData class

#include "main.hpp" 
#include "test.hpp"

#include "mesh.hpp"
#include "data.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;


  //----------------------------------------------------------------------
  unit_init(0,1);
  //----------------------------------------------------------------------

  FieldDescr * field_descr = new FieldDescr;

  int i1 = field_descr->insert_permanent("f1");
  int i2 = field_descr->insert_permanent("f2");
  int i3 = field_descr->insert_permanent("f3");
  int i4 = field_descr->insert_permanent("f4");
  int i5 = field_descr->insert_permanent("f5");

  int it1 = field_descr->insert_temporary("t1");
  int it2 = field_descr->insert_temporary("t2");
  int it3 = field_descr->insert_temporary("t3");

  // set precision

  field_descr->set_precision(i1, precision_single);
  field_descr->set_precision(i2, precision_double);
  field_descr->set_precision(i3, precision_double);
  field_descr->set_precision(i4, precision_double);
  field_descr->set_precision(i5, precision_quadruple);
  field_descr->set_precision(it1, precision_single);
  field_descr->set_precision(it2, precision_double);
  field_descr->set_precision(it3, precision_double);

  unit_class ("Cello");

  unit_func ("32-bit floating-point");
  unit_assert (sizeof(float)       == 4);
  unit_func ("64-bit floating-point");
  unit_assert (sizeof(double)      == 8);
  unit_func ("128-bit floating-point");
  unit_assert (sizeof(long double) == 16);

  bool is_quad_supported = sizeof(long double) == 16;

  // set ghosts

  int g1[3] = {1,1,1};
  int g2[3] = {2,2,2};
  int g3[3] = {1,0,2};
  int g4[3] = {0,0,0};
  int g5[3] = {3,3,3};

  field_descr->set_ghost_depth(i1, g1[0],g1[1],g1[2]);
  field_descr->set_ghost_depth(i2, g2[0],g2[1],g2[2]);
  field_descr->set_ghost_depth(i3, g3[0],g3[1],g3[2]);
  field_descr->set_ghost_depth(i4, g4[0],g4[1],g4[2]);
  field_descr->set_ghost_depth(i5, g5[0],g5[1],g5[2]);
  field_descr->set_ghost_depth(it1, g1[0],g1[1],g1[2]);
  field_descr->set_ghost_depth(it2, g2[0],g2[1],g2[2]);
  field_descr->set_ghost_depth(it3, g3[0],g3[1],g3[2]);

  // set centering

  field_descr->set_centering(i2, 1, 0, 0);
  field_descr->set_centering(i3, 0, 1, 0);
  field_descr->set_centering(i4, 0, 0, 1);

  field_descr->set_centering(it1, 1, 0, 0);
  field_descr->set_centering(it2, 0, 1, 0);
  field_descr->set_centering(it3, 0, 0, 1);

  int nx,ny,nz;
  nx=4; ny=8; nz=6;
  double xpm,ypm,zpm;
  xpm = -1.0;  ypm = -2.0, zpm = -3.0;
  double xpp,ypp,zpp;
  xpp =  1.0;  ypp =  2.0, zpp =  3.0;

  FieldData * field_data = new FieldData(field_descr, nx,ny,nz);
  //----------------------------------------------------------------------

  unit_class("FieldData");

  //----------------------------------------------------------------------

  unit_func("size");

  int size[3];

  field_data->size(size+0,size+1,size+2);

  unit_assert(size[0]==nx && size[1]==ny && size[2]==nz);

  //----------------------------------------------------------------------
  // allocate / deallocate
  //----------------------------------------------------------------------

  unit_func("permanent_allocated");
  unit_assert( ! field_data->permanent_allocated());

  // Allocate

  unit_func("allocate_permanent");

  field_data->allocate_permanent(field_descr,false);

  size_t array_size_without_ghosts = field_data->permanent_size();

  unit_assert(field_data->permanent() != 0);
  unit_assert(field_data->permanent_allocated());
  unit_assert(field_data->permanent_size() > 0);

  // Reallocate

  unit_func("reallocate_permanent");

  field_data->reallocate_permanent(field_descr,true);

  size_t array_size_with_ghosts = field_data->permanent_size();

  unit_assert(field_data->permanent() != 0);
  unit_assert(field_data->permanent_allocated());
  unit_assert(field_data->permanent_size() > 0);
  unit_assert(array_size_with_ghosts > array_size_without_ghosts);

  // Deallocate

  unit_func("deallocate_permanent");

  field_data->deallocate_permanent();

  unit_assert(field_data->permanent() == 0);
  unit_assert( ! field_data->permanent_allocated());
  unit_assert(field_data->permanent_size() == 0);

  // Allocate

  unit_func("allocate_permanent");

  field_data->allocate_permanent(field_descr,true);

  unit_assert(field_data->permanent() != 0);
  unit_assert(field_data->permanent_allocated());
  unit_assert(field_data->permanent_size() == array_size_with_ghosts);

  //----------------------------------------------------------------------
  field_data->reallocate_permanent(field_descr,false);

  float       *v1,*u1;
  double      *v2,*u2;
  double      *v3,*u3;
  double      *v4,*u4;
  long double *v5,*u5;
  float       *vt1;
  double      *vt2;
  double      *vt3;

  unit_func("values");  // without ghosts

  v1 = (float *)       field_data->values(field_descr,i1);
  v2 = (double *)      field_data->values(field_descr,i2);
  v3 = (double *)      field_data->values(field_descr,i3);
  v4 = (double *)      field_data->values(field_descr,i4);
  v5 = (long double *) field_data->values(field_descr,i5);

  unit_func("temporary values");

  vt1 = (float *)   field_data->values(field_descr,it1);
  vt2 = (double *)  field_data->values(field_descr,it2);
  vt3 = (double *)  field_data->values(field_descr,it3);
  unit_assert (vt1 == 0);
  unit_assert (vt2 == 0);
  unit_assert (vt3 == 0);
  field_data->allocate_temporary(field_descr,it1);
  field_data->allocate_temporary(field_descr,it2);
  field_data->allocate_temporary(field_descr,it3);
  vt1 = (float *)   field_data->values(field_descr,it1);
  vt2 = (double *)  field_data->values(field_descr,it2);
  vt3 = (double *)  field_data->values(field_descr,it3);
  unit_assert (vt1 != 0);
  unit_assert (vt2 != 0);
  unit_assert (vt3 != 0);
  field_data->deallocate_temporary(field_descr,it1);
  field_data->deallocate_temporary(field_descr,it2);
  field_data->deallocate_temporary(field_descr,it3);
  vt1 = (float *)   field_data->values(field_descr,it1);
  vt2 = (double *)  field_data->values(field_descr,it2);
  vt3 = (double *)  field_data->values(field_descr,it3);
  unit_assert (vt1 == 0);
  unit_assert (vt2 == 0);
  unit_assert (vt3 == 0);

  // field sizes without ghosts

  int nu1 = (nx)   * (ny)   * (nz);
  int nu2 = (nx+1) * (ny)   * (nz);
  int nu3 = (nx)   * (ny+1) * (nz);
  int nu4 = (nx)   * (ny)   * (nz+1);

  // field sizes with ghosts

  int nv1 = (nx + 2*g1[0])   * (ny + 2*g1[1])   * (nz + 2*g1[2]);
  int nv2 = (nx + 2*g2[0]+1) * (ny + 2*g2[1])   * (nz + 2*g2[2]);
  int nv3 = (nx + 2*g3[0])   * (ny + 2*g3[1]+1) * (nz + 2*g3[2]);
  int nv4 = (nx + 2*g4[0])   * (ny + 2*g4[1])   * (nz + 2*g4[2]+1);
  
  size_t nb1 = (char *) v2 - (char *) v1;
  size_t nb2 = (char *) v3 - (char *) v2;
  size_t nb3 = (char *) v4 - (char *) v3;
  size_t nb4 = (char *) v5 - (char *) v4;

  unit_assert (nb1 == sizeof (float) * nu1);
  printf("nb2,nu2 = %lu %lu",nb2,sizeof(double)*nu2);
  unit_assert (nb2 == sizeof (double)* nu2);
  printf("nb3,nu3 = %lu %lu",nb3,sizeof(double)*nu3);
  unit_assert (nb3 == sizeof (double)* nu3);
  unit_assert (nb4 == sizeof (double)* nu4);

  //----------------------------------------------------------------------

  unit_func("unknowns");  // without ghosts

  typedef float type_f1;
  u1 = (type_f1 *)      field_data->unknowns(field_descr,i1);
  u2 = (double *)       field_data->unknowns(field_descr,i2);
  u3 = (double *)       field_data->unknowns(field_descr,i3);
  u4 = (double *)       field_data->unknowns(field_descr,i4);
  u5 = (long double *)  field_data->unknowns(field_descr,i5);

  unit_assert(u1 != 0);
  unit_assert(u2 != 0);
  unit_assert(u3 != 0);
  unit_assert(u4 != 0);
  unit_assert(u5 != 0);

  nb1 = (char *)u2 - (char *)u1;
  nb2 = (char *)u3 - (char *)u2;
  nb3 = (char *)u4 - (char *)u3;
  nb4 = (char *)u5-(char *)u4;

  unit_assert (nb1 == sizeof (float) * nu1);
  unit_assert (nb2 == sizeof (double)* nu2);
  unit_assert (nb3 == sizeof (double)* nu3);
  unit_assert (nb4 == sizeof (double)* nu4);

  // unknown field
  unit_assert (field_data->unknowns (field_descr,1000) == NULL);
  unit_assert (field_data->unknowns (field_descr,"not a field") == NULL);


  //----------------------------------------------------------------------

  // with ghosts

  field_data->reallocate_permanent(field_descr,true);

  v1 =  (float *)      field_data->values(field_descr,i1);
  v2 = (double *)      field_data->values(field_descr,i2);
  v3 = (double *)      field_data->values(field_descr,i3);
  v4 = (double *)      field_data->values(field_descr,i4);
  v5 = (long double *) field_data->values(field_descr,i5);
  
  unit_assert(field_data->values(field_descr,i1) == 
	      field_data->values(field_descr,"f1"));
  unit_assert(field_data->values(field_descr,i2) == 
	      field_data->values(field_descr,"f2"));
  unit_assert(field_data->values(field_descr,i3) == 
	      field_data->values(field_descr,"f3"));
  unit_assert(field_data->values(field_descr,i4) == 
	      field_data->values(field_descr,"f4"));
  unit_assert(field_data->values(field_descr,i5) == 
	      field_data->values(field_descr,"f5"));

  unit_assert(v1 != 0);
  unit_assert(v2 != 0);
  unit_assert(v3 != 0);
  unit_assert(v4 != 0);
  unit_assert(v5 != 0);

  nb1 = (char *)v2 - (char *)v1;
  nb2 = (char *)v3 - (char *)v2;
  nb3 = (char *)v4 - (char *)v3;
  nb4 = (char *)v5 - (char *)v4;

  unit_assert (nb1 == sizeof (float)  * nv1);
  unit_assert (nb2 == sizeof (double) * nv2);
  unit_assert (nb3 == sizeof (double) * nv3);
  unit_assert (nb4 == sizeof (double) * nv4);

  // unknown field
  unit_assert (field_data->values (field_descr,1000) == NULL);
  unit_assert (field_data->values (field_descr,"not a field") == NULL);

  unit_func("unknowns");  // with ghosts

  u1 = (float *)       field_data->unknowns(field_descr,i1);
  u2 = (double *)      field_data->unknowns(field_descr,i2);
  u3 = (double *)      field_data->unknowns(field_descr,i3);
  u4 = (double *)      field_data->unknowns(field_descr,i4);
  u5 = (long double *) field_data->unknowns(field_descr,i5);

  unit_assert(u1 != 0);
  unit_assert(u2 != 0);
  unit_assert(u3 != 0);
  unit_assert(u4 != 0);
  unit_assert(u5 != 0);

  nb1 = (char *)u2 - (char *)u1;
  nb2 = (char *)u3 - (char *)u2;
  nb3 = (char *)u4 - (char *)u3;
  nb4 = (char *)u5 - (char *)u4;

  // a,b fields  u unknowns  v values  g ghosts
  // bu - au = (bv + bg) - (av + ag) 
  //         = (bv + (bu-bv)) - (av + (au-av))
  //         = (bv - av) + (bu-bv) - (au-av)


  int n1[3] = { nx+2*g1[0],   ny+2*g1[1],   nz+2*g1[2]};
  int n2[3] = { nx+2*g2[0]+1, ny+2*g2[1],   nz+2*g2[2]};
  int n3[3] = { nx+2*g3[0],   ny+2*g3[1]+1, nz+2*g3[2]};
  int n4[3] = { nx+2*g4[0],   ny+2*g4[1],   nz+2*g4[2]+1};
  int n5[3] = { nx+2*g5[0],   ny+2*g5[1],   nz+2*g5[2]};

  size_t ng1 = sizeof (float)      * (g1[0] + n1[0] *(g1[1] + n1[1]*g1[2]));
  size_t ng2 = sizeof (double)     * (g2[0] + n2[0] *(g2[1] + n2[1]*g2[2]));
  size_t ng3 = sizeof (double)     * (g3[0] + n3[0] *(g3[1] + n3[1]*g3[2]));
  size_t ng4 = sizeof (double)     * (g4[0] + n4[0] *(g4[1] + n4[1]*g4[2]));
  size_t ng5 = sizeof (long double)* (g5[0] + n5[0] *(g5[1] + n5[1]*g5[2]));

  unit_assert (nb1 == sizeof (float) * nv1 + ng2 - ng1);
  unit_assert (nb2 == sizeof (double)* nv2 + ng3 - ng2);
  unit_assert (nb3 == sizeof (double)* nv3 + ng4 - ng3);
  if (is_quad_supported) unit_assert (nb4 == sizeof (double)* nv4 + ng5 - ng4);

  bool passed;

  //--------------------------------------------------
  // Fill values interior then test against unknowns

  for (int iz=0; iz<n1[2]; iz++) {
    for (int iy=0; iy<n1[1]; iy++) {
      for (int ix=0; ix<n1[0]; ix++) {
	int i = ix + n1[0]*(iy + n1[1]*iz);
	v1[i] = 
	  (g1[0] <= ix && ix < n1[0]-g1[0] &&
	   g1[1] <= iy && iy < n1[1]-g1[1] &&
	   g1[2] <= iz && iz < n1[2]-g1[2]) ? 10 : -10;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n1[0]*(iy + n1[1]*iz);
	if (u1[i] != 10) passed = false;
      }
    }
  }

  unit_assert(passed);

  //--------------------------------------------------

  for (int iz=0; iz<n2[2]; iz++) {
    for (int iy=0; iy<n2[1]; iy++) {
      for (int ix=0; ix<n2[0]; ix++) {
	int i = ix + n2[0]*(iy + n2[1]*iz);
	v2[i] = 
	  (g2[0] <= ix && ix < n2[0]-g2[0] &&
	   g2[1] <= iy && iy < n2[1]-g2[1] &&
	   g2[2] <= iz && iz < n2[2]-g2[2]) ? 20 : -20;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n2[0]*(iy + n2[1]*iz);
	if (u2[i] != 20) passed = false;
      }
    }
  }

  unit_assert(passed);

  //--------------------------------------------------

  for (int iz=0; iz<n3[2]; iz++) {
    for (int iy=0; iy<n3[1]; iy++) {
      for (int ix=0; ix<n3[0]; ix++) {
	int i = ix + n3[0]*(iy + n3[1]*iz);
	v3[i] = 
	  (g3[0] <= ix && ix < n3[0]-g3[0] &&
	   g3[1] <= iy && iy < n3[1]-g3[1] &&
	   g3[2] <= iz && iz < n3[2]-g3[2]) ? 30 : -30;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n3[0]*(iy + n3[1]*iz);
	if (u3[i] != 30) passed = false;
      }
    }
  }

  unit_assert(passed);

  //--------------------------------------------------

  for (int iz=0; iz<n4[2]; iz++) {
    for (int iy=0; iy<n4[1]; iy++) {
      for (int ix=0; ix<n4[0]; ix++) {
	int i = ix + n4[0]*(iy + n4[1]*iz);
	v4[i] = 
	  (g4[0] <= ix && ix < n4[0]-g4[0] &&
	   g4[1] <= iy && iy < n4[1]-g4[1] &&
	   g4[2] <= iz && iz < n4[2]-g4[2]) ? 40 : -40;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n4[0]*(iy + n4[1]*iz);
	if (u4[i] != 40) passed = false;
      }
    }
  }

  unit_assert(passed);

  //--------------------------------------------------


  for (int iz=0; iz<n5[2]; iz++) {
    for (int iy=0; iy<n5[1]; iy++) {
      for (int ix=0; ix<n5[0]; ix++) {
	int i = ix + n5[0]*(iy + n5[1]*iz);
	v5[i] = 
	  (g5[0] <= ix && ix < n5[0]-g5[0] &&
	   g5[1] <= iy && iy < n5[1]-g5[1] &&
	   g5[2] <= iz && iz < n5[2]-g5[2]) ? 50 : -50;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n5[0]*(iy + n5[1]*iz);
	if (u5[i] != 50) passed = false;
      }
    }
  }

  if (is_quad_supported) unit_assert(passed);

  //--------------------------------------------------

  // test that first value after f1 is f2 ghost, etc.

  unit_assert(*((double *)(v1+n1[0]*n1[1]*n1[2]))     == -20);
  unit_assert(*((double *)(v2+n2[0]*n2[1]*n2[2]))     == -30);
  unit_assert(*((double *)(v3+n3[0]*n3[1]*n3[2]))     ==  40);  // g4 has no ghosts
  unit_assert(*((long double *)(v4+n4[0]*n4[1]*n4[2])) == -50);

  unit_assert(*((float *) (v2-1)) == -10);
  unit_assert(*((double *) (v3-1)) == -20);
  unit_assert(*((double *) (v4-1)) == -30);
  unit_assert(*((double *) (v5-1)) ==  40); // g4 has no ghosts

  //----------------------------------------------------------------------
  unit_func("cell_width");

  double hx=0,hy=0,hz=0;

  field_data->cell_width(xpm,xpp,&hx,ypm,ypp,&hy,zpm,zpp,&hz);

  unit_assert(fabs(hx-2.0/nx) < 1e-6);
  unit_assert(fabs(hy-4.0/ny) < 1e-6);
  unit_assert(fabs(hz-6.0/nz) < 1e-6);

	
  //----------------------------------------------------------------------
  unit_func("clear");

  field_data->clear(field_descr,2.0);

  unit_assert(2.0 == v1[0]);
  unit_assert(2.0 == v5[n5[0]*n5[1]*n5[2]-1]);

  field_data->clear(field_descr,3.0, 1 );
  
  unit_assert(2.0 == v1[n1[0]*n1[1]*n1[2]-1]);
  unit_assert(3.0 == v2[0] );
  unit_assert(3.0 == v2[n2[0]*n2[1]*n2[2]-1]);
  unit_assert(2.0 == v3[0] );
	
  field_data->clear(field_descr,4.0, 2, 3 );

  unit_assert(3.0 == v2[n2[0]*n2[1]*n2[2]-1]);
  unit_assert(4.0 == v3[0] );
  unit_assert(4.0 == v3[n3[0]*n3[1]*n3[2]-1]);
  unit_assert(4.0 == v4[0] );
  unit_assert(4.0 == v4[n4[0]*n4[1]*n4[2]-1]);
  unit_assert(2.0 == v5[0] );

  //----------------------------------------------------------------------
  unit_func("reallocate_ghosts");

  field_data->reallocate_permanent(field_descr,false);

  v1 = (float *) field_data->values(field_descr,i1);
  v2 = (double *) field_data->values(field_descr,i2);
  v3 = (double *) field_data->values(field_descr,i3);
  v4 = (double *) field_data->values(field_descr,i4);
  v5 = (long double *) field_data->values(field_descr,i5);

  unit_assert(3.0 == v2[(nx+1)*ny*nz-1]);
  unit_assert(4.0 == v3[0] );
  unit_assert(4.0 == v3[nx*(ny+1)*nz-1]);
  unit_assert(4.0 == v4[0] );
  unit_assert(4.0 == v4[nx*ny*(nz+1)-1]);
  unit_assert(2.0 == v5[0] );
  
  //----------------------------------------------------------------------
  unit_finalize();
  //----------------------------------------------------------------------

  exit_();
}
PARALLEL_MAIN_END
