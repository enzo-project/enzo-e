// See LICENSE_CELLO file for license and copyright information

/// @file     test_Field.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the Field class

#include "main.hpp" 
#include "test.hpp"

#include "mesh.hpp"
#include "data.hpp"

struct field_info_type {
  int field_density;
  int field_velocity_x;
  int field_velocity_y;
  int field_velocity_z;
  int field_total_energy;

  int group_density;
  int group_vector;

  int gx, gy, gz;
  int cx, cy, cz;
};

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;


  //----------------------------------------------------------------------
  unit_init(0,1);
  //----------------------------------------------------------------------

  //----------------------------------------------------------------------
  // FieldData tests from test_FieldData
  //----------------------------------------------------------------------
  {

    int nx,ny,nz;
    nx=4; ny=5; nz=6;
    FieldDescr * field_descr = new FieldDescr;
    FieldData * field_data = new FieldData(field_descr, nx,ny,nz);

    Field * field = new Field(field_descr,field_data);
  
    int i1 = field->insert_permanent("f1");
    int i2 = field->insert_permanent("f2");
    int i3 = field->insert_permanent("f3");
    int i4 = field->insert_permanent("f4");
    int i5 = field->insert_permanent("f5");

    // set precision

    field->set_precision(i1, precision_single);
    field->set_precision(i2, precision_double);
    field->set_precision(i3, precision_double);
    field->set_precision(i4, precision_double);
    field->set_precision(i5, precision_quadruple);

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

    field->set_ghosts(i1, g1[0],g1[1],g1[2]);
    field->set_ghosts(i2, g2[0],g2[1],g2[2]);
    field->set_ghosts(i3, g3[0],g3[1],g3[2]);
    field->set_ghosts(i4, g4[0],g4[1],g4[2]);
    field->set_ghosts(i5, g5[0],g5[1],g5[2]);

    // set centering

    field->set_centering(i2, 1, 0, 0);
    field->set_centering(i3, 0, 1, 0);
    field->set_centering(i4, 0, 0, 1);

    double xpm,ypm,zpm;
    xpm = -1.0;  ypm = -2.0, zpm = -3.0;
    double xpp,ypp,zpp;
    xpp =  1.0;  ypp =  2.0, zpp =  3.0;


    //----------------------------------------------------------------------

    unit_class("FieldData");

    //----------------------------------------------------------------------

    unit_func("size");

    int size[3];

    field->size(size+0,size+1,size+2);

    unit_assert(size[0]==nx && size[1]==ny && size[2]==nz);

    //----------------------------------------------------------------------
    // allocate / deallocate
    //----------------------------------------------------------------------

    unit_func("array_allocated");
    unit_assert( ! field->permanent_allocated());

    // Allocate

    unit_func("allocate_permanent");

    field->allocate_permanent(false);

    size_t array_size_without_ghosts = field->permanent_size();

    unit_assert(field->permanent() != 0);
    unit_assert(field->permanent_allocated());
    unit_assert(field->permanent_size() > 0);

    // Reallocate

    unit_func("reallocate_permanent");

    field->reallocate_permanent(true);

    size_t array_size_with_ghosts = field->permanent_size();

    unit_assert(field->permanent() != 0);
    unit_assert(field->permanent_allocated());
    unit_assert(field->permanent_size() > 0);
    unit_assert(array_size_with_ghosts > array_size_without_ghosts);

    // Deallocate

    unit_func("deallocate_permanent");

    field->deallocate_permanent();

    unit_assert(field->permanent() == 0);
    unit_assert( ! field->permanent_allocated());
    unit_assert(field->permanent_size() == 0);

    // Allocate

    unit_func("allocate_permanent");

    field->allocate_permanent(true);

    unit_assert(field->permanent() != 0);
    unit_assert(field->permanent_allocated());
    unit_assert(field->permanent_size() == array_size_with_ghosts);

  
    //----------------------------------------------------------------------
    field->reallocate_permanent(false);


    float       *v1,*u1;
    double      *v2,*u2;
    double      *v3,*u3;
    double      *v4,*u4;
    long double *v5,*u5;

    unit_func("values");  // without ghosts
  
    v1 = (float *)       field->values(i1);
    v2 = (double *)      field->values(i2);
    v3 = (double *)      field->values(i3);
    v4 = (double *)      field->values(i4);
    v5 = (long double *) field->values(i5);
  
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
    TRACE2("nb2,nu2 = %d %d",nb2,sizeof(double)*nu2);
    unit_assert (nb2 == sizeof (double)* nu2);
    TRACE2("nb3,nu3 = %d %d",nb3,sizeof(double)*nu3);
    unit_assert (nb3 == sizeof (double)* nu3);
    unit_assert (nb4 == sizeof (double)* nu4);

    //----------------------------------------------------------------------

    unit_func("unknowns");  // without ghosts

    typedef float type_f1;
    u1 = (type_f1 *)      field->unknowns(i1);
    u2 = (double *)       field->unknowns(i2);
    u3 = (double *)       field->unknowns(i3);
    u4 = (double *)       field->unknowns(i4);
    u5 = (long double *)  field->unknowns(i5);

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
    unit_assert (field->unknowns (1000) == NULL);
    unit_assert (field->unknowns ("not a field") == NULL);


    //----------------------------------------------------------------------

    // with ghosts

    field->reallocate_permanent(true);

    v1 =  (float *)      field->values(i1);
    v2 = (double *)      field->values(i2);
    v3 = (double *)      field->values(i3);
    v4 = (double *)      field->values(i4);
    v5 = (long double *) field->values(i5);
  
    unit_assert(field->values(i1) == 
		field->values("f1"));
    unit_assert(field->values(i2) == 
		field->values("f2"));
    unit_assert(field->values(i3) == 
		field->values("f3"));
    unit_assert(field->values(i4) == 
		field->values("f4"));
    unit_assert(field->values(i5) == 
		field->values("f5"));

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
    unit_assert (field->values (1000) == NULL);
    unit_assert (field->values ("not a field") == NULL);

    unit_func("unknowns");  // with ghosts

    u1 = (float *)       field->unknowns(i1);
    u2 = (double *)      field->unknowns(i2);
    u3 = (double *)      field->unknowns(i3);
    u4 = (double *)      field->unknowns(i4);
    u5 = (long double *) field->unknowns(i5);

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

    field->cell_width(xpm,xpp,&hx,ypm,ypp,&hy,zpm,zpp,&hz);

    unit_assert(fabs(hx-2.0/nx) < 1e-6);
    unit_assert(fabs(hy-4.0/ny) < 1e-6);
    unit_assert(fabs(hz-6.0/nz) < 1e-6);

	
    //----------------------------------------------------------------------
    unit_func("clear");

    field->clear(2.0);

    unit_assert(2.0 == v1[0]);
    unit_assert(2.0 == v5[n5[0]*n5[1]*n5[2]-1]);

    field->clear(3.0, 1 );
  
    unit_assert(2.0 == v1[n1[0]*n1[1]*n1[2]-1]);
    unit_assert(3.0 == v2[0] );
    unit_assert(3.0 == v2[n2[0]*n2[1]*n2[2]-1]);
    unit_assert(2.0 == v3[0] );
	
    field->clear(4.0, 2, 3 );

    unit_assert(3.0 == v2[n2[0]*n2[1]*n2[2]-1]);
    unit_assert(4.0 == v3[0] );
    unit_assert(4.0 == v3[n3[0]*n3[1]*n3[2]-1]);
    unit_assert(4.0 == v4[0] );
    unit_assert(4.0 == v4[n4[0]*n4[1]*n4[2]-1]);
    unit_assert(2.0 == v5[0] );

    //----------------------------------------------------------------------
    unit_func("reallocate_ghosts");

    field->reallocate_permanent(false);

    v1    = 
      (float *) field->values(i1);
    v2 = 
				 (double *) field->values(i2);
    v3 = 
      (double *) field->values(i3);
    v4 = 
      (double *) field->values(i4);
    v5 =
      (long double *) field->values(i5);

    unit_assert(3.0 == v2[(nx+1)*ny*nz-1]);
    unit_assert(4.0 == v3[0] );
    unit_assert(4.0 == v3[nx*(ny+1)*nz-1]);
    unit_assert(4.0 == v4[0] );
    unit_assert(4.0 == v4[nx*ny*(nz+1)-1]);
    unit_assert(2.0 == v5[0] );

    unit_func("delete");
    delete field->field_descr();
    delete field->field_data();
    delete field;
  }

  //----------------------------------------------------------------------
  // FieldDescr Tests from test_FieldDescr
  //----------------------------------------------------------------------

  {
    unit_class("FieldDescr");

    struct field_info_type info;

    unit_func ("FieldDescr");
    FieldDescr * field_descr = new FieldDescr;
    FieldData * field_data = 0;

    Field * field = new Field(field_descr,field_data);

    unit_assert(field_descr != 0);
    printf ("sizeof(FieldDescr) = %lud\n",sizeof(FieldDescr));

    // Fields

    unit_func("insert_permanent");
    unit_assert(field->field_count()==0);
    field->insert_permanent("density");
    unit_assert(field->field_count()==1);
    field->insert_permanent("velocity_x");
    unit_assert(field->field_count()==2);
    field->insert_permanent("velocity_y");
    unit_assert(field->field_count()==3);
    field->insert_permanent("velocity_z");
    unit_assert(field->field_count()==4);
    field->insert_permanent("total_energy");
    unit_assert(field->field_count()==5);
    field->insert_permanent("total_energy");
    unit_assert(field->field_count()==5);

    unit_func("field_count");
    unit_assert(field->field_count()==5);

    unit_func("field_id");

    info.field_density      = field->field_id("density");
    info.field_velocity_x   = field->field_id("velocity_x");
    info.field_velocity_y   = field->field_id("velocity_y");
    info.field_velocity_z   = field->field_id("velocity_z");
    info.field_total_energy = field->field_id("total_energy");

    unit_assert(field->field_id("density")      == info.field_density);
    unit_assert(field->field_id("velocity_x")   == info.field_velocity_x);
    unit_assert(field->field_id("velocity_y")   == info.field_velocity_y);
    unit_assert(field->field_id("velocity_z")   == info.field_velocity_z);
    unit_assert(field->field_id("total_energy") == info.field_total_energy);

    unit_func("is_field");

    unit_assert(field->is_field("density"));
    unit_assert(! field->is_field("not_a_field"));

    unit_func("field_name");

    unit_assert(field->field_name(info.field_density)      == "density");
    unit_assert(field->field_name(info.field_velocity_x)   == "velocity_x");
    unit_assert(field->field_name(info.field_velocity_y)   == "velocity_y");
    unit_assert(field->field_name(info.field_velocity_z)   == "velocity_z");
    unit_assert(field->field_name(info.field_total_energy) == "total_energy");

    //----------------------------------------------------------------------
    // Global attributes
    //----------------------------------------------------------------------

    // (set and reset in case test value is a default)

    field->set_alignment(8);
    field->set_padding(64);
    field->set_courant(0.5);
  
    unit_func("alignment");
    unit_assert(field->alignment() == 8);
    unit_func("padding");
    unit_assert(field->padding() == 64);
    unit_func("courant");
    unit_assert(field->courant() == 0.5);

    field->set_alignment(4);
    field->set_padding(32);
    field->set_courant(0.75);
  
    unit_func("alignment");
    unit_assert(field->alignment() == 4);
    unit_func("padding");
    unit_assert(field->padding() == 32);
    unit_func("courant");
    unit_assert(field->courant() == 0.75);
  
    //----------------------------------------------------------------------
    // Field attributes
    //----------------------------------------------------------------------

    // Precision

    unit_func("precision");

    field->set_precision(info.field_density,    precision_single);
    field->set_precision(info.field_velocity_x, precision_double);
    field->set_precision(info.field_velocity_y, precision_double);
    field->set_precision(info.field_velocity_z, precision_double);

    unit_assert(field->precision(info.field_density)      == precision_single);
    unit_assert(field->precision(info.field_velocity_x)   == precision_double);
    unit_assert(field->precision(info.field_velocity_y)   == precision_double);
    unit_assert(field->precision(info.field_velocity_z)   == precision_double);
    unit_assert(field->precision(info.field_total_energy) == default_precision);
  

    unit_func("bytes_per_element");
    unit_assert(field->bytes_per_element(info.field_density)==4);

    // Centering

    unit_func("centering");

    field->set_centering(info.field_velocity_x, 1, 0, 0);
    field->set_centering(info.field_velocity_y, 0, 1, 0);
    field->set_centering(info.field_velocity_z, 0, 0, 1);


    field->centering(info.field_density, &info.cx, &info.cy, &info.cz);
    unit_assert(info.cx == 0 && info.cy == 0 && info.cz == 0);

    field->centering(info.field_velocity_x, &info.cx, &info.cy, &info.cz);
    unit_assert(info.cx == 1 && info.cy == 0 && info.cz == 0);

    field->centering(info.field_velocity_y, &info.cx, &info.cy, &info.cz);
    unit_assert(info.cx == 0 && info.cy == 1 && info.cz == 0);

    field->centering(info.field_velocity_z, &info.cx, &info.cy, &info.cz);
    unit_assert(info.cx == 0 && info.cy == 0 && info.cz == 1);
  
    // Ghost zone depth

    unit_func("ghosts");

    field->set_ghosts(info.field_density, 3, 3, 3);
    field->set_ghosts(info.field_velocity_x, 1, 0, 0);
    field->set_ghosts(info.field_velocity_y, 0, 1, 0);
    field->set_ghosts(info.field_velocity_z, 0, 0, 1);

    field->ghosts(info.field_density, &info.gx, &info.gy, &info.gz);
    unit_assert(info.gx==3 && info.gy==3 && info.gz==3);
    field->ghosts(info.field_velocity_x, &info.gx, &info.gy, &info.gz);
    unit_assert(info.gx==1 && info.gy==0 && info.gz==0);
    field->ghosts(info.field_velocity_y, &info.gx, &info.gy, &info.gz);
    unit_assert(info.gx==0 && info.gy==1 && info.gz==0);
    field->ghosts(info.field_velocity_z, &info.gx, &info.gy, &info.gz);
    unit_assert(info.gx==0 && info.gy==0 && info.gz==1);


  }
  //----------------------------------------------------------------------
  unit_finalize();
  //----------------------------------------------------------------------

  exit_();
}
PARALLEL_MAIN_END
