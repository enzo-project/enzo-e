// See LICENSE_CELLO file for license and copyright information

/// @file     test_Field.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the Field class

#include "main.hpp" 
#include "test.hpp"

#include "mesh.hpp"
#include "view.hpp"
#include "data.hpp"

/// helper function that actually checks consistency between the view and the
/// field data
template<class T, class FieldT>
bool check_view_props_(CelloView<T,3> view, FieldT& field, int id,
                       bool coarse, bool include_ghosts,
                       int history = 0){
  int nx, ny, nz;
  T* ptr;

  if (coarse){
    ASSERT("check_view_props_",
           "require history==0 and include_ghosts==true for coarse fields",
           (history == 0) && include_ghosts);
    ptr = reinterpret_cast<T*>(field.coarse_values(id));
    field.coarse_dimensions(id, &nx, &ny, &nz);
  } else {
    field.dimensions(id, &nx, &ny, &nz);
    if (include_ghosts & (!field.ghosts_allocated())){
      ERROR("check_view_props_", "something's wrong ghosts must be allocated");
    } else if (include_ghosts) {
      ptr = reinterpret_cast<T*>(field.values(id, history));
    } else {
      ptr = reinterpret_cast<T*>(field.unknowns(id, history));
      int gx, gy, gz;
      field.ghost_depth(id, &gx, &gy, &gz);
      if (nx > 1){ nx -= 2*gx; }
      if (ny > 1){ ny -= 2*gy; }
      if (nz > 1){ nz -= 2*gz; }
    }
  }

  bool matching_shape = ((view.shape(0) == nz) &
                         (view.shape(1) == ny) &
                         (view.shape(2) == nx));
  bool passed = (view.data() == ptr) & (ptr != nullptr) & matching_shape;
  return passed;
}

/// Tests the basic properties of the views of a field.
///
/// This namely just checks the view's shape and the address of the first
/// element. This function tries to exhaust many permutation of options to
/// create a view.
template<class T>
void test_view_(Field& field, int id, std::string name = "",
                bool coarse = false, int history = 0){

  bool passed = true;

  if (coarse){
    ASSERT("test_view",
           "require history==0 and name==\"\" for coarse fields",
           (history == 0) && (name == ""));
    { // first, handle non-const Field stuff
      CelloView<T, 3> view = field.coarse_view<T>(id);
      passed &= check_view_props_(view, field, id, true, true, history);
    }
    { // now handle const Field stuff
      const Field& const_field = field;
      CelloView<const T, 3> view = const_field.coarse_view<T>(id);
      passed &= check_view_props_(view, const_field, id, true, true, history);
    }

  } else {

    auto run_test = [&](ghost_choice gchoice, bool includes_ghost)
      {
        // first handle non-const Field stuff
        // load from id
        {
          CelloView<T, 3> view = field.view<T>(id, gchoice, history);
          passed &= check_view_props_(view, field, id, false, includes_ghost,
                                      history);
        }
        // load from field name (if it was provided)
        if (name != ""){
          CelloView<T, 3> view = field.view<T>(name, gchoice, history);
          passed &= check_view_props_(view, field, id, false, includes_ghost,
                                      history);
        }

        // now, handle const Field
        {
          const Field& const_field = field;
          CelloView<const T, 3> view =
            const_field.view<T>(id, gchoice, history);
          passed &= check_view_props_(view, const_field, id, false,
                                      includes_ghost, history);
        }
      };

    // first, test ghost_choice::exclude
    run_test(ghost_choice::exclude, false);

    // now, test ghost_choice::permit
    run_test(ghost_choice::permit, field.ghosts_allocated());

    // finally, if ghost zones are allocated, check ghost_choice::include
    if (field.ghosts_allocated()){
      run_test(ghost_choice::include, true);
    }
  }
  unit_assert(passed);
}

struct field_info_type {
  int field_density;
  int field_velocity_x;
  int field_velocity_y;
  int field_velocity_z;
  int field_total_energy;


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
    nx=4; ny=8; nz=6;
    FieldDescr * field_descr = new FieldDescr;
    FieldData * field_data = new FieldData(field_descr, nx,ny,nz);

    Field field(field_descr,field_data);
  
    int i1 = field.insert_permanent("f1");
    int i2 = field.insert_permanent("f2");
    int i3 = field.insert_permanent("f3");
    int i4 = field.insert_permanent("f4");
    int i5 = field.insert_permanent("f5");

    int j1 = field.insert_temporary();
    int j2 = field.insert_temporary();
    int j3 = field.insert_temporary("j3");

    // set precision

    field.set_precision(i1, precision_single);
    field.set_precision(i2, precision_double);
    field.set_precision(i3, precision_double);
    field.set_precision(i4, precision_double);
    field.set_precision(i5, precision_quadruple);

    field.set_precision(j1,precision_single);
    field.set_precision(j2,precision_double);
    field.set_precision(j3,precision_quadruple);

    unit_class("Field");
    unit_func("num_permanent");
    unit_assert(field.num_permanent() == 5);
    unit_func("num_temporary");
    unit_assert(field.num_temporary() == 3);
    
    unit_class ("Cello");

    unit_func ("32-bit floating-point");
    unit_assert (sizeof(float)       == 4);
    unit_func ("64-bit floating-point");
    unit_assert (sizeof(double)      == 8);
    unit_func ("128-bit floating-point");
    unit_assert (sizeof(long double) == 16);

    bool have_quad = sizeof(long double) == 16;

    // set ghosts

    int g1[3] = {1,1,1};
    int g2[3] = {2,2,2};
    int g3[3] = {1,0,2};
    int g4[3] = {0,0,0};
    int g5[3] = {3,3,3};

    field.set_ghost_depth(i1, g1[0],g1[1],g1[2]);
    field.set_ghost_depth(i2, g2[0],g2[1],g2[2]);
    field.set_ghost_depth(i3, g3[0],g3[1],g3[2]);
    field.set_ghost_depth(i4, g4[0],g4[1],g4[2]);
    field.set_ghost_depth(i5, g5[0],g5[1],g5[2]);

    // set centering

    field.set_centering(i2, 1, 0, 0);
    field.set_centering(i3, 0, 1, 0);
    field.set_centering(i4, 0, 0, 1);

    double xpm,ypm,zpm;
    xpm = -1.0;  ypm = -2.0, zpm = -3.0;
    double xpp,ypp,zpp;
    xpp =  1.0;  ypp =  2.0, zpp =  3.0;

    // set history

    unit_func("set_history");
    
    field.set_history (3);

    unit_assert(field.num_history() == 3);

    //----------------------------------------------------------------------

    unit_class("FieldData");

    //----------------------------------------------------------------------

    unit_func("size");

    int size[3];

    field.size(size+0,size+1,size+2);

    unit_assert(size[0]==nx && size[1]==ny && size[2]==nz);

    //----------------------------------------------------------------------
    // allocate / deallocate
    //----------------------------------------------------------------------

    unit_func("array_allocated");
    unit_assert( ! field.permanent_allocated());

    // Allocate

    unit_func("allocate_permanent");

    field.allocate_permanent(false);

    size_t array_size_without_ghosts = field.permanent_size();

    unit_assert(field.permanent() != 0);
    unit_assert(field.permanent_allocated());
    unit_assert(field.permanent_size() > 0);

    // Reallocate

    unit_func("reallocate_permanent");

    field.reallocate_permanent(true);

    size_t array_size_with_ghosts = field.permanent_size();

    unit_assert(field.permanent() != 0);
    unit_assert(field.permanent_allocated());
    unit_assert(field.permanent_size() > 0);
    unit_assert(array_size_with_ghosts > array_size_without_ghosts);

    // Deallocate

    unit_func("deallocate_permanent");

    field.deallocate_permanent();

    unit_assert(field.permanent() == 0);
    unit_assert( ! field.permanent_allocated());
    unit_assert(field.permanent_size() == 0);

    // Allocate

    unit_func("allocate_permanent");

    field.allocate_permanent(true);

    unit_assert(field.permanent() != 0);
    unit_assert(field.permanent_allocated());
    unit_assert(field.permanent_size() == array_size_with_ghosts);

  
    //----------------------------------------------------------------------
    field.reallocate_permanent(false);


    float       *v1,*u1;
    double      *v2,*u2;
    double      *v3,*u3;
    double      *v4,*u4;
    long double *v5,*u5;

    unit_func("values");  // without ghosts
  
    v1 = (float *)       field.values(i1);
    v2 = (double *)      field.values(i2);
    v3 = (double *)      field.values(i3);
    v4 = (double *)      field.values(i4);
    v5 = (long double *) field.values(i5);
  
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
    u1 = (type_f1 *)      field.unknowns(i1);
    u2 = (double *)       field.unknowns(i2);
    u3 = (double *)       field.unknowns(i3);
    u4 = (double *)       field.unknowns(i4);
    u5 = (long double *)  field.unknowns(i5);

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
    unit_assert (field.unknowns (1000) == NULL);
    unit_assert (field.unknowns ("not a field") == NULL);


    //----------------------------------------------------------------------

    // with ghosts

    field.reallocate_permanent(true);

    v1 =  (float *)      field.values(i1);
    v2 = (double *)      field.values(i2);
    v3 = (double *)      field.values(i3);
    v4 = (double *)      field.values(i4);
    v5 = (long double *) field.values(i5);
  
    unit_assert(field.values(i1) == 
		field.values("f1"));
    unit_assert(field.values(i2) == 
		field.values("f2"));
    unit_assert(field.values(i3) == 
		field.values("f3"));
    unit_assert(field.values(i4) == 
		field.values("f4"));
    unit_assert(field.values(i5) == 
		field.values("f5"));

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
    unit_assert (field.values (1000) == NULL);
    unit_assert (field.values ("not a field") == NULL);

    unit_func("unknowns");  // with ghosts

    u1 = (float *)       field.unknowns(i1);
    u2 = (double *)      field.unknowns(i2);
    u3 = (double *)      field.unknowns(i3);
    u4 = (double *)      field.unknowns(i4);
    u5 = (long double *) field.unknowns(i5);

    unit_assert(u1 != 0);
    unit_assert(u2 != 0);
    unit_assert(u3 != 0);
    unit_assert(u4 != 0);
    unit_assert(u5 != 0);

    nb1 = (char *)u2 - (char *)u1;
    nb2 = (char *)u3 - (char *)u2;
    nb3 = (char *)u4 - (char *)u3;
    nb4 = (char *)u5 - (char *)u4;

    // test value arrays differ from unknown arrays by ghost offset
    // a,b fields  u unknowns  v values  g ghosts
    // bu - au = (bv + bg) - (av + ag) 
    //         = (bv + (bu-bv)) - (av + (au-av))
    //         = (bv - av) + (bu-bv) - (au-av)

    int m1[3] = { nx+2*g1[0],   ny+2*g1[1],   nz+2*g1[2]};
    int m2[3] = { nx+2*g2[0]+1, ny+2*g2[1],   nz+2*g2[2]};
    int m3[3] = { nx+2*g3[0],   ny+2*g3[1]+1, nz+2*g3[2]};
    int m4[3] = { nx+2*g4[0],   ny+2*g4[1],   nz+2*g4[2]+1};
    int m5[3] = { nx+2*g5[0],   ny+2*g5[1],   nz+2*g5[2]};

    size_t ng1 = sizeof (float)      * (g1[0] + m1[0] *(g1[1] + m1[1]*g1[2]));
    size_t ng2 = sizeof (double)     * (g2[0] + m2[0] *(g2[1] + m2[1]*g2[2]));
    size_t ng3 = sizeof (double)     * (g3[0] + m3[0] *(g3[1] + m3[1]*g3[2]));
    size_t ng4 = sizeof (double)     * (g4[0] + m4[0] *(g4[1] + m4[1]*g4[2]));
    size_t ng5 = sizeof (long double)* (g5[0] + m5[0] *(g5[1] + m5[1]*g5[2]));

    unit_assert (nb1 == sizeof (float) * nv1 + ng2 - ng1);
    unit_assert (nb2 == sizeof (double)* nv2 + ng3 - ng2);
    unit_assert (nb3 == sizeof (double)* nv3 + ng4 - ng3);
    if (have_quad) unit_assert (nb4 == sizeof (double)* nv4 + ng5 - ng4);

    //--------------------------------------------------
    // Fill values interior then test against unknowns

    for (int iz=0; iz<m1[2]; iz++) {
      for (int iy=0; iy<m1[1]; iy++) {
	for (int ix=0; ix<m1[0]; ix++) {
	  int i = ix + m1[0]*(iy + m1[1]*iz);
	  v1[i] = 
	    (g1[0] <= ix && ix < m1[0]-g1[0] &&
	     g1[1] <= iy && iy < m1[1]-g1[1] &&
	     g1[2] <= iz && iz < m1[2]-g1[2]) ? 10 : -10;

	}
      }
    }

    bool passed = true;

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int i = ix + m1[0]*(iy + m1[1]*iz);
	  if (u1[i] != 10) passed = false;
	}
      }
    }

    unit_assert(passed);

    //--------------------------------------------------

    for (int iz=0; iz<m2[2]; iz++) {
      for (int iy=0; iy<m2[1]; iy++) {
	for (int ix=0; ix<m2[0]; ix++) {
	  int i = ix + m2[0]*(iy + m2[1]*iz);
	  v2[i] = 
	    (g2[0] <= ix && ix < m2[0]-g2[0] &&
	     g2[1] <= iy && iy < m2[1]-g2[1] &&
	     g2[2] <= iz && iz < m2[2]-g2[2]) ? 20 : -20;

	}
      }
    }

    passed = true;

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int i = ix + m2[0]*(iy + m2[1]*iz);
	  if (u2[i] != 20) passed = false;
	}
      }
    }

    unit_assert(passed);

    //--------------------------------------------------

    for (int iz=0; iz<m3[2]; iz++) {
      for (int iy=0; iy<m3[1]; iy++) {
	for (int ix=0; ix<m3[0]; ix++) {
	  int i = ix + m3[0]*(iy + m3[1]*iz);
	  v3[i] = 
	    (g3[0] <= ix && ix < m3[0]-g3[0] &&
	     g3[1] <= iy && iy < m3[1]-g3[1] &&
	     g3[2] <= iz && iz < m3[2]-g3[2]) ? 30 : -30;

	}
      }
    }

    passed = true;

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int i = ix + m3[0]*(iy + m3[1]*iz);
	  if (u3[i] != 30) passed = false;
	}
      }
    }

    unit_assert(passed);

    //--------------------------------------------------

    for (int iz=0; iz<m4[2]; iz++) {
      for (int iy=0; iy<m4[1]; iy++) {
	for (int ix=0; ix<m4[0]; ix++) {
	  int i = ix + m4[0]*(iy + m4[1]*iz);
	  v4[i] = 
	    (g4[0] <= ix && ix < m4[0]-g4[0] &&
	     g4[1] <= iy && iy < m4[1]-g4[1] &&
	     g4[2] <= iz && iz < m4[2]-g4[2]) ? 40 : -40;

	}
      }
    }

    passed = true;

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int i = ix + m4[0]*(iy + m4[1]*iz);
	  if (u4[i] != 40) passed = false;
	}
      }
    }

    unit_assert(passed);

    //--------------------------------------------------


    for (int iz=0; iz<m5[2]; iz++) {
      for (int iy=0; iy<m5[1]; iy++) {
	for (int ix=0; ix<m5[0]; ix++) {
	  int i = ix + m5[0]*(iy + m5[1]*iz);
	  v5[i] = 
	    (g5[0] <= ix && ix < m5[0]-g5[0] &&
	     g5[1] <= iy && iy < m5[1]-g5[1] &&
	     g5[2] <= iz && iz < m5[2]-g5[2]) ? 50 : -50;

	}
      }
    }

    passed = true;

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int i = ix + m5[0]*(iy + m5[1]*iz);
	  if (u5[i] != 50) passed = false;
	}
      }
    }

    if (have_quad) unit_assert(passed);

    //--------------------------------------------------

    // test that first value after f1 is f2 ghost, etc.

    const int M1 = m1[0]*m1[1]*m1[2];
    const int M2 = m2[0]*m2[1]*m2[2];
    const int M3 = m3[0]*m3[1]*m3[2];
    const int M4 = m4[0]*m4[1]*m4[2];
    const int M5 = m5[0]*m5[1]*m5[2];

    unit_assert(*((double *)(v1+M1))     == -20);
    unit_assert(*((double *)(v2+M2))     == -30);
    unit_assert(*((double *)(v3+M3))     ==  40);
    unit_assert(*((long double *)(v4+M4)) == -50);

    unit_assert(*((float *) (v2-1)) == -10);
    unit_assert(*((double *) (v3-1)) == -20);
    unit_assert(*((double *) (v4-1)) == -30);
    unit_assert(*((double *) (v5-1)) ==  40);

    //----------------------------------------------------------------------
    unit_func("cell_width");

    double hx=0,hy=0,hz=0;

    field.cell_width(xpm,xpp,&hx,ypm,ypp,&hy,zpm,zpp,&hz);

    unit_assert(fabs(hx-2.0/nx) < 1e-6);
    unit_assert(fabs(hy-4.0/ny) < 1e-6);
    unit_assert(fabs(hz-6.0/nz) < 1e-6);

	
    //----------------------------------------------------------------------
    unit_func("clear");

    field.clear(2.0);

    unit_assert(2.0 == v1[0]);
    unit_assert(2.0 == v5[M5-1]);

    field.clear(3.0, 1 );
  
    unit_assert(2.0 == v1[M1-1]);
    unit_assert(3.0 == v2[0] );
    unit_assert(3.0 == v2[M2-1]);
    unit_assert(2.0 == v3[0] );
	
    field.clear(4.0, 2, 3 );

    unit_assert(3.0 == v2[M2-1]);
    unit_assert(4.0 == v3[0] );
    unit_assert(4.0 == v3[M3-1]);
    unit_assert(4.0 == v4[0] );
    unit_assert(4.0 == v4[M4-1]);
    unit_assert(2.0 == v5[0] );

    //----------------------------------------------------------------------
    unit_func("reallocate_ghosts");

    field.reallocate_permanent(false);

    v1 =       (float *) field.values(i1);
    v2 =      (double *) field.values(i2);
    v3 =      (double *) field.values(i3);
    v4 =      (double *) field.values(i4);
    v5 = (long double *) field.values(i5);

    unit_assert(3.0 == v2[(nx+1)*ny*nz-1]);
    unit_assert(4.0 == v3[0] );
    unit_assert(4.0 == v3[nx*(ny+1)*nz-1]);
    unit_assert(4.0 == v4[0] );
    unit_assert(4.0 == v4[nx*ny*(nz+1)-1]);
    unit_assert(2.0 == v5[0] );

    // Test temporary fields
    
    float * t1;
    double * t2;
    long double * t3;

    unit_func("insert_temporary");

    t1 = (float *) field.values(j1);
    t2 = (double *) field.values(j2);
    t3 = (long double *) field.values("j3");

    unit_assert (t1 == NULL);
    unit_assert (t2 == NULL);
    unit_assert (t3 == NULL);
    
    unit_func("allocate_temporary");

    field.allocate_temporary(j1);    
    field.allocate_temporary(j2);    
    field.allocate_temporary("j3");

    t1 = (float *)       field.values(j1);
    t2 = (double *)      field.values(j2);
    t3 = (long double *) field.values("j3");

    t1[0] = -1.00;
    t2[0] = -2.00;
    t3[0] = -3.00;

    unit_assert ((t1[0] == -1.00));
    unit_assert ((t2[0] == -2.00));
    unit_assert ((t3[0] == -3.00));

    // Test coarse fields

    // unit_assert (cv1 == NULL);
    // unit_assert (cv2 == NULL);
    // unit_assert (cv3 == NULL);
    // unit_assert (cv4 == NULL);
    // unit_assert (cv5 == NULL);
    // unit_assert (ct1 == NULL);
    // unit_assert (ct2 == NULL);
    // unit_assert (ct3 == NULL);
    
    // unit_func("allocate_coarse");

    field.allocate_coarse();

    float       *cv1 = (float *) field.coarse_values(i1);
    double      *cv2 = (double *) field.coarse_values(i2);
    double      *cv3 = (double *) field.coarse_values(i3);
    double      *cv4 = (double *) field.coarse_values(i4);
    long double *cv5 = (long double *) field.coarse_values(i5);
    float       *ct1 = (float *) field.coarse_values(j1);
    double      *ct2 = (double *) field.coarse_values(j2);
    double      *ct3 = (double *) field.coarse_values(j3);

    cv1[0] = -1.00;
    cv2[0] = -2.00;
    cv3[0] = -3.00;
    cv4[0] = -4.00;
    cv5[0] = -5.00;
    ct1[0] = 1.00;
    ct2[0] = 2.00;
    ct3[0] = 3.00;

    unit_assert ((cv1[0] == -1.00));
    unit_assert ((cv2[0] == -2.00));
    unit_assert ((cv3[0] == -3.00));
    unit_assert ((cv4[0] == -4.00));
    unit_assert ((cv5[0] == -5.00));
    unit_assert ((ct1[0] == 1.00));
    unit_assert ((ct2[0] == 2.00));
    unit_assert ((ct3[0] == 3.00));

    // ---------------------------------------------------------------------
    unit_func("view");

    for (int i = 0; i < 2; i++){
      for (int history_ind = 0; history_ind < 3; history_ind++){
        field.reallocate_permanent(i == 0);

        // permanent fields
        test_view_<float>(field, i1, "f1", false, history_ind);
        test_view_<double>(field, i2, "f2", false, history_ind);
        test_view_<double>(field, i3, "f3", false, history_ind);
        test_view_<double>(field, i4, "f4", false, history_ind);
        test_view_<long double>(field, i5, "f5", false, history_ind);

        // temporary fields
        test_view_<float>(field, j1, "", false, history_ind);
        test_view_<double>(field, j2, "", false, history_ind);
        test_view_<long double>(field, j3, "j3", false, history_ind);

      }
    }

    // ---------------------------------------------------------------------
    unit_func("coarse_view");

    for (int i = 0; i < 2; i++){
      field.reallocate_permanent(i == 0);

      // permanent fields
      test_view_<float>(field, i1, "", true, 0);
      test_view_<double>(field, i2, "", true, 0);
      test_view_<double>(field, i3, "", true, 0);
      test_view_<double>(field, i4, "", true, 0);
      test_view_<long double>(field, i5, "", true, 0);

      test_view_<float>(field, j1, "", false, 0);
      test_view_<double>(field, j2, "", false, 0);
      test_view_<long double>(field, j3, "", false, 0);
    }

    //--------------------------------------------------

    // History

#define HIST_INIT(ip,ih,i) (17.0*(ip+3) + 7.0*(ih+2) - 13.0*(i+7))
    
    field.reallocate_permanent(true);

    v1 =       (float *) field.values(i1);
    v2 =      (double *) field.values(i2);
    v3 =      (double *) field.values(i3);
    v4 =      (double *) field.values(i4);
    v5 = (long double *) field.values(i5);

    // Initialize and save what will be history 3

    int ih = 3;
    
    for (int i=0; i<M1; i++) {  v1[i] = HIST_INIT(1,ih,i); }
    for (int i=0; i<M2; i++) {  v2[i] = HIST_INIT(2,ih,i); }
    for (int i=0; i<M3; i++) {  v3[i] = HIST_INIT(3,ih,i); }
    for (int i=0; i<M4; i++) {  v4[i] = HIST_INIT(4,ih,i); }
    for (int i=0; i<M5; i++) {  v5[i] = HIST_INIT(5,ih,i); }

    field.save_history(5.0-ih);

    // Initialize and save what will be history 2

    ih = 2;

    for (int i=0; i<M1; i++) {  v1[i] = HIST_INIT(1,ih,i); }
    for (int i=0; i<M2; i++) {  v2[i] = HIST_INIT(2,ih,i); }
    for (int i=0; i<M3; i++) {  v3[i] = HIST_INIT(3,ih,i); }
    for (int i=0; i<M4; i++) {  v4[i] = HIST_INIT(4,ih,i); }
    for (int i=0; i<M5; i++) {  v5[i] = HIST_INIT(5,ih,i); }

    field.save_history(5.0-ih);

    // Initialize and save what will be history 1

    ih = 1;

    for (int i=0; i<M1; i++) {  v1[i] = HIST_INIT(1,ih,i); }
    for (int i=0; i<M2; i++) {  v2[i] = HIST_INIT(2,ih,i); }
    for (int i=0; i<M3; i++) {  v3[i] = HIST_INIT(3,ih,i); }
    for (int i=0; i<M4; i++) {  v4[i] = HIST_INIT(4,ih,i); }
    for (int i=0; i<M5; i++) {  v5[i] = HIST_INIT(5,ih,i); }

    field.save_history(5.0-ih);

    // Initialize current values (index_history = 0)
    
    ih = 0;

    for (int i=0; i<M1; i++) {  v1[i] = HIST_INIT(1,ih,i); }
    for (int i=0; i<M2; i++) {  v2[i] = HIST_INIT(2,ih,i); }
    for (int i=0; i<M3; i++) {  v3[i] = HIST_INIT(3,ih,i); }
    for (int i=0; i<M4; i++) {  v4[i] = HIST_INIT(4,ih,i); }
    for (int i=0; i<M5; i++) {  v5[i] = HIST_INIT(5,ih,i); }

    float       *v1h1,*v1h2,*v1h3;
    double      *v2h1,*v2h2,*v2h3;
    double      *v3h1,*v3h2,*v3h3;
    double      *v4h1,*v4h2,*v4h3;
    long double *v5h1,*v5h2,*v5h3;

    v1h1 = (float *) field.values(i1,1);
    v1h2 = (float *) field.values(i1,2);
    v1h3 = (float *) field.values(i1,3);
    v2h1 = (double *) field.values(i2,1);
    v2h2 = (double *) field.values(i2,2);
    v2h3 = (double *) field.values(i2,3);
    v3h1 = (double *) field.values(i3,1);
    v3h2 = (double *) field.values(i3,2);
    v3h3 = (double *) field.values(i3,3);
    v4h1 = (double *) field.values(i4,1);
    v4h2 = (double *) field.values(i4,2);
    v4h3 = (double *) field.values(i4,3);
    v5h1 = (long double *) field.values(i5,1);
    v5h2 = (long double *) field.values(i5,2);
    v5h3 = (long double *) field.values(i5,3);

    // test current values
    ih = 0;
    unit_func ("history[0]");
    passed = true;
    for (int i=0; i<M1; i++) { passed &= (v1[i] == HIST_INIT(1,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M2; i++) { passed &= (v2[i] == HIST_INIT(2,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M3; i++) { passed &= (v3[i] == HIST_INIT(3,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M4; i++) { passed &= (v4[i] == HIST_INIT(4,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M5; i++) { passed &= (v5[i] == HIST_INIT(5,ih,i));  }
    unit_assert (passed);

    // note--no time saved for history = 0

    ih = 1;
    unit_func ("history[1]");
    passed = true;
    for (int i=0; i<M1; i++) { passed &= (v1h1[i] == HIST_INIT(1,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M2; i++) { passed &= (v2h1[i] == HIST_INIT(2,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M3; i++) { passed &= (v3h1[i] == HIST_INIT(3,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M4; i++) { passed &= (v4h1[i] == HIST_INIT(4,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M5; i++) { passed &= (v5h1[i] == HIST_INIT(5,ih,i));  }
    unit_assert (passed);

    unit_func ("history_time[1]");
    unit_assert (5.0-ih == field.history_time(ih-1));

    ih = 2;
    unit_func ("history[2]");
    passed = true;
    for (int i=0; i<M1; i++) { passed &= (v1h2[i] == HIST_INIT(1,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M2; i++) { passed &= (v2h2[i] == HIST_INIT(2,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M3; i++) { passed &= (v3h2[i] == HIST_INIT(3,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M4; i++) { passed &= (v4h2[i] == HIST_INIT(4,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M5; i++) { passed &= (v5h2[i] == HIST_INIT(5,ih,i));  }
    unit_assert (passed);

    unit_func ("history_time[2]");
    unit_assert (5.0-ih == field.history_time(ih-1));

    ih = 3;
    unit_func ("history[3]");
    passed = true;
    for (int i=0; i<M1; i++) { passed &= (v1h3[i] == HIST_INIT(1,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M2; i++) { passed &= (v2h3[i] == HIST_INIT(2,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M3; i++) { passed &= (v3h3[i] == HIST_INIT(3,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M4; i++) { passed &= (v4h3[i] == HIST_INIT(4,ih,i));  }
    unit_assert (passed);
    passed = true;
    for (int i=0; i<M5; i++) { passed &= (v5h3[i] == HIST_INIT(5,ih,i));  }
    unit_assert (passed);

    unit_func ("history_time[3]");
    unit_assert (5.0-ih == field.history_time(ih-1));

    //--------------------------------------------------
    unit_func ("units_scale_cgs");

    double s1 = 2.0;
    double s2 = 4.0;
    double s3 = 8.0;
    double s4 = 16.0;
    double s5 = 32.0;

    // initialize fields in code units
    for (int i=0; i<M1; i++) {  v1[i] = HIST_INIT(1,0,i); }
    for (int i=0; i<M2; i++) {  v2[i] = HIST_INIT(2,0,i); }
    for (int i=0; i<M3; i++) {  v3[i] = HIST_INIT(3,0,i); }
    for (int i=0; i<M4; i++) {  v4[i] = HIST_INIT(4,0,i); }
    for (int i=0; i<M5; i++) {  v5[i] = HIST_INIT(5,0,i); }

    // scale fields to cgs given scaling

    field.units_scale_cgs (i1,s1);
    field.units_scale_cgs (i2,s2);
    field.units_scale_cgs (i3,s3);
    field.units_scale_cgs (i4,s4);
    field.units_scale_cgs (i5,s5);

    // test fields are scaled
    passed = true;
    for (int i=0; i<M1; i++) { passed &= (v1[i] == s1*HIST_INIT(1,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M2; i++) { passed &= (v2[i] == s2*HIST_INIT(2,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M3; i++) { passed &= (v3[i] == s3*HIST_INIT(3,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M4; i++) { passed &= (v4[i] == s4*HIST_INIT(4,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M5; i++) { passed &= (v5[i] == s5*HIST_INIT(5,0,i)); }
    unit_assert(passed);

    // try to rescale fields again

    field.units_scale_cgs (i1,s1);
    field.units_scale_cgs (i2,s2);
    field.units_scale_cgs (i3,s3);
    field.units_scale_cgs (i4,s4);
    field.units_scale_cgs (i5,s5);

    // check that they were not scaled twice

    passed = true;
    for (int i=0; i<M1; i++) { passed &= (v1[i] == s1*HIST_INIT(1,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M2; i++) { passed &= (v2[i] == s2*HIST_INIT(2,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M3; i++) { passed &= (v3[i] == s3*HIST_INIT(3,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M4; i++) { passed &= (v4[i] == s4*HIST_INIT(4,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M5; i++) { passed &= (v5[i] == s5*HIST_INIT(5,0,i)); }
    unit_assert(passed);

    // unscale fields back to code units

    unit_func ("units_scale_code");

    field.units_scale_code (i1,s1);
    field.units_scale_code (i2,s2);
    field.units_scale_code (i3,s3);
    field.units_scale_code (i4,s4);
    field.units_scale_code (i5,s5);

    // test that fields are in original code units
    
    passed = true;
    for (int i=0; i<M1; i++) { passed &= (v1[i] == HIST_INIT(1,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M2; i++) { passed &= (v2[i] == HIST_INIT(2,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M3; i++) { passed &= (v3[i] == HIST_INIT(3,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M4; i++) { passed &= (v4[i] == HIST_INIT(4,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M5; i++) { passed &= (v5[i] == HIST_INIT(5,0,i)); }
    unit_assert(passed);

    // try to re-scale fields to code units

    field.units_scale_code (i1,s1);
    field.units_scale_code (i2,s2);
    field.units_scale_code (i3,s3);
    field.units_scale_code (i4,s4);
    field.units_scale_code (i5,s5);

    // check that fields were not unscaled twice

    passed = true;
    for (int i=0; i<M1; i++) { passed &= (v1[i] == HIST_INIT(1,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M2; i++) { passed &= (v2[i] == HIST_INIT(2,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M3; i++) { passed &= (v3[i] == HIST_INIT(3,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M4; i++) { passed &= (v4[i] == HIST_INIT(4,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M5; i++) { passed &= (v5[i] == HIST_INIT(5,0,i)); }
    unit_assert(passed);

    // convert to cgs using units  s

    unit_func ("units_scale_cgs");

    field.units_scale_cgs (i1,s1);
    field.units_scale_cgs (i2,s2);
    field.units_scale_cgs (i3,s3);
    field.units_scale_cgs (i4,s4);
    field.units_scale_cgs (i5,s5);

    // reconvert to cgs using updated units r

    double r1 = 0.25;
    double r2 = 0.125;
    double r3 = 64.0;
    double r4 = 32.0;
    double r5 = 16.0;

    field.units_scale_cgs (i1,r1);
    field.units_scale_cgs (i2,r2);
    field.units_scale_cgs (i3,r3);
    field.units_scale_cgs (i4,r4);
    field.units_scale_cgs (i5,r5);

    // test that fields were unscaled by s first before being re-scaled by r

    passed = true;
    for (int i=0; i<M1; i++) { passed &= (v1[i] == r1*HIST_INIT(1,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M2; i++) { passed &= (v2[i] == r2*HIST_INIT(2,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M3; i++) { passed &= (v3[i] == r3*HIST_INIT(3,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M4; i++) { passed &= (v4[i] == r4*HIST_INIT(4,0,i)); }
    unit_assert(passed);
    passed = true;
    for (int i=0; i<M5; i++) { passed &= (v5[i] == r5*HIST_INIT(5,0,i)); }
    unit_assert(passed);

    
    //--------------------------------------------------
    
    unit_func("deallocate_temporary");

    field.deallocate_temporary(j1);    
    field.deallocate_temporary(j2);    
    field.deallocate_temporary("j3");

    t1 = (float *) field.values(j1);
    t2 = (double *) field.values(j2);
    t3 = (long double *) field.values("j3");

    unit_assert (t1 == NULL);
    unit_assert (t2 == NULL);
    unit_assert (t3 == NULL);
    
    unit_func("delete");
    delete field.field_descr();
    delete field.field_data();
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

    Field field(field_descr,field_data);

    unit_assert(field_descr != 0);
    printf ("sizeof(FieldDescr) = %lud\n",sizeof(FieldDescr));

    // Fields

    unit_func("insert_permanent");
    unit_assert(field.field_count()==0);
    field.insert_permanent("density");
    unit_assert(field.field_count()==1);
    field.insert_permanent("velocity_x");
    unit_assert(field.field_count()==2);
    field.insert_permanent("velocity_y");
    unit_assert(field.field_count()==3);
    field.insert_permanent("velocity_z");
    unit_assert(field.field_count()==4);
    field.insert_permanent("total_energy");
    unit_assert(field.field_count()==5);
    field.insert_permanent("total_energy");
    unit_assert(field.field_count()==5);


    unit_func("field_count");
    unit_assert(field.field_count()==5);

    unit_func("field_id");

    info.field_density      = field.field_id("density");
    info.field_velocity_x   = field.field_id("velocity_x");
    info.field_velocity_y   = field.field_id("velocity_y");
    info.field_velocity_z   = field.field_id("velocity_z");
    info.field_total_energy = field.field_id("total_energy");

    unit_assert(field.field_id("density")      == info.field_density);
    unit_assert(field.field_id("velocity_x")   == info.field_velocity_x);
    unit_assert(field.field_id("velocity_y")   == info.field_velocity_y);
    unit_assert(field.field_id("velocity_z")   == info.field_velocity_z);
    unit_assert(field.field_id("total_energy") == info.field_total_energy);

    unit_func("is_field");

    unit_assert(field.is_field("density"));
    unit_assert(! field.is_field("not_a_field"));

    unit_func("field_name");

    unit_assert(field.field_name(info.field_density)      == "density");
    unit_assert(field.field_name(info.field_velocity_x)   == "velocity_x");
    unit_assert(field.field_name(info.field_velocity_y)   == "velocity_y");
    unit_assert(field.field_name(info.field_velocity_z)   == "velocity_z");
    unit_assert(field.field_name(info.field_total_energy) == "total_energy");

    //----------------------------------------------------------------------
    // Global attributes
    //----------------------------------------------------------------------

    // (set and reset in case test value is a default)

    field.set_alignment(8);
    field.set_padding(64);
  
    unit_func("alignment");
    unit_assert(field.alignment() == 8);
    unit_func("padding");
    unit_assert(field.padding() == 64);

    field.set_alignment(4);
    field.set_padding(32);
  
    unit_func("alignment");
    unit_assert(field.alignment() == 4);
    unit_func("padding");
    unit_assert(field.padding() == 32);
  
    //----------------------------------------------------------------------
    // Field attributes
    //----------------------------------------------------------------------

    // Precision

    unit_func("precision");

    field.set_precision(info.field_density,    precision_single);
    field.set_precision(info.field_velocity_x, precision_double);
    field.set_precision(info.field_velocity_y, precision_double);
    field.set_precision(info.field_velocity_z, precision_double);

    unit_assert(field.precision(info.field_density)      == precision_single);
    unit_assert(field.precision(info.field_velocity_x)   == precision_double);
    unit_assert(field.precision(info.field_velocity_y)   == precision_double);
    unit_assert(field.precision(info.field_velocity_z)   == precision_double);
    unit_assert(field.precision(info.field_total_energy) == default_precision);
  
    unit_func("bytes_per_element");
    unit_assert(field.bytes_per_element(info.field_density)==4);
    unit_assert(field.bytes_per_element(info.field_velocity_x)==8);
    unit_assert(field.bytes_per_element(info.field_velocity_y)==8);
    unit_assert(field.bytes_per_element(info.field_velocity_z)==8);

    // Centering

    unit_func("centering");

    field.set_centering(info.field_velocity_x, 1, 0, 0);
    field.set_centering(info.field_velocity_y, 0, 1, 0);
    field.set_centering(info.field_velocity_z, 0, 0, 1);


    field.centering(info.field_density, &info.cx, &info.cy, &info.cz);
    unit_assert(info.cx == 0 && info.cy == 0 && info.cz == 0);

    field.centering(info.field_velocity_x, &info.cx, &info.cy, &info.cz);
    unit_assert(info.cx == 1 && info.cy == 0 && info.cz == 0);

    field.centering(info.field_velocity_y, &info.cx, &info.cy, &info.cz);
    unit_assert(info.cx == 0 && info.cy == 1 && info.cz == 0);

    field.centering(info.field_velocity_z, &info.cx, &info.cy, &info.cz);
    unit_assert(info.cx == 0 && info.cy == 0 && info.cz == 1);
  
    // Ghost zone depth

    unit_func("ghosts");

    field.set_ghost_depth(info.field_density, 3, 3, 3);
    field.set_ghost_depth(info.field_velocity_x, 1, 0, 0);
    field.set_ghost_depth(info.field_velocity_y, 0, 1, 0);
    field.set_ghost_depth(info.field_velocity_z, 0, 0, 1);

    field.ghost_depth(info.field_density, &info.gx, &info.gy, &info.gz);
    unit_assert(info.gx==3 && info.gy==3 && info.gz==3);
    field.ghost_depth(info.field_velocity_x, &info.gx, &info.gy, &info.gz);
    unit_assert(info.gx==1 && info.gy==0 && info.gz==0);
    field.ghost_depth(info.field_velocity_y, &info.gx, &info.gy, &info.gz);
    unit_assert(info.gx==0 && info.gy==1 && info.gz==0);
    field.ghost_depth(info.field_velocity_z, &info.gx, &info.gy, &info.gz);
    unit_assert(info.gx==0 && info.gy==0 && info.gz==1);

  }
  //----------------------------------------------------------------------
  unit_finalize();
  //----------------------------------------------------------------------

  exit_();
}
PARALLEL_MAIN_END
