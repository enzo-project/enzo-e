// See LICENSE_CELLO file for license and copyright information

/// @file     method_InitialDefault.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-02-16
/// @todo     Add image mask for field initialization
/// @todo     Parse initialization parameters once per Simulation rather than once per Block
/// @brief    Implementation of the InitialDefault class
///
/// Detailed description of file method_InitialDefault.cpp

#include "cello.hpp"

#include "method.hpp"

//----------------------------------------------------------------------

InitialDefault::InitialDefault(Parameters * parameters) throw ()
  : Initial (),
    parameters_(parameters)
{
}

//----------------------------------------------------------------------

void InitialDefault::compute (const FieldDescr * field_descr,
			      Block * block) throw()
{
  // Initialize Fields according to parameters

  //--------------------------------------------------
  parameters_->group_set(0,"Initial");
  //--------------------------------------------------

  FieldBlock *       field_block = block->field_block();

  double *value=0, *vdeflt=0, *x=0, *y=0, *z=0, *t=0;
  bool * region=0, *rdeflt=0;
  int n, nx=0,ny=0,nz=0;

  for (int index_field = 0;
       index_field < field_descr->field_count();
       index_field++) {
    
    std::string field_name = field_descr->field_name(index_field);

    //--------------------------------------------------
    parameters_->group_set(1,field_name);
    //--------------------------------------------------

    // If Initial:<field_name>:value is a list, try parsing it

    // parameter: Initial : <field> : value

    parameter_enum parameter_type = parameters_->type("value");

    if (parameter_type == parameter_float) {

      field_block->clear(field_descr,parameters_->value_float("value",0.0), 
			 index_field);

    } else if (parameter_type == parameter_list) {

      // Check parameter length

      int list_length = parameters_->list_length("value");

      ASSERT1("InitialDefault::compute",
	     "Length of list parameter Initial:%s:value must be odd",
	     field_name.c_str(),
	     (list_length % 2) == 1);

      // Allocate arrays if needed
      if (value == NULL) {
	allocate_xyzt_(block,index_field,
		       &nx,&ny,&nz,
		       &value, &vdeflt,
		       &region,&rdeflt,
		       &x,&y,&z,&t);
      }

      // Evaluate last non-conditional equation in list

      n = nx*ny*nz;

      evaluate_float_ (field_block, list_length-1, field_name, 
			n, value,vdeflt,x,y,z,t);

      for (int i=0; i<n; i++) region[i] = true;

      copy_values_ (field_descr,field_block,value, region,index_field,nx,ny,nz);

      // Evaluate conditional equations in list

      for (int index_list=0; index_list < list_length-1; index_list+=2) {

	evaluate_float_ (field_block, index_list, field_name, 
			  n, value,vdeflt,x,y,z,t);

	evaluate_logical_ (field_block, index_list + 1, field_name, 
			  n, region,rdeflt,x,y,z,t);

	copy_values_ (field_descr,field_block,value, region,index_field,nx,ny,nz);

      }

    } else if (parameter_type == parameter_unknown) {
      WARNING1("InitialDefault::compute",  
	       "Uninitialized field %s",
	       field_name.c_str());
      
    } else {
      ERROR2("InitialDefault::compute",
	     "Illegal parameter type %s when initializing field %s",
	     parameter_type_name[parameter_type],field_name.c_str());
    }
  }
  // Deallocate arrays if needed
  if (value != NULL) {
    delete [] value;
    delete [] vdeflt;
    delete [] region;
    delete [] rdeflt;
    delete [] x;
    delete [] y;
    delete [] z;
    delete [] t;
  }
}

//======================================================================

void InitialDefault::allocate_xyzt_
(
 Block * block,
 int index_field,
 int * nx, int * ny, int * nz,
 double ** value, double ** vdeflt,
 bool   ** region,bool   ** rdeflt,
 double ** x, double ** y, double ** z, double ** t
 ) throw()
{

  FieldBlock * field_block = block->field_block();

  // Get field size

  field_block->size(nx,ny,nz);

  int n = (*nx)*(*ny)*(*nz);

  // Allocate double arrays

  (*value)  = new double [n];
  (*vdeflt) = new double [n];
  (*region) = new bool [n];
  (*rdeflt) = new bool [n];
  (*x)      = new double [n];
  (*y)      = new double [n];
  (*z)      = new double [n];
  (*t)      = new double [n];

  double xm, xp, ym, yp, zm, zp;

  block->lower(&xm,&ym,&zm);
  block->upper(&xp,&yp,&zp);

  double hx,hy,hz;
  field_block->cell_width(block,&hx,&hy,&hz);

  // Initialize arrays
  for (int iz=0; iz<(*nz); iz++) {
    for (int iy=0; iy<(*ny); iy++) {
      for (int ix=0; ix<(*nx); ix++) {
	int i=ix + (*nx)*(iy + (*ny)*iz);
	(*value)[i] = 0.0;
	(*vdeflt)[i]  = 0.0;
	(*region)[i] = false;
	(*rdeflt)[i]  = false;
	(*x)[i] = xm + ix*hx;
	(*y)[i] = ym + iy*hy;
	(*z)[i] = zm + iz*hz;
	(*t)[i] = 0.0;
	//	printf ("%d %d %d %g %g %g\n",ix,iy,iz,(*x)[i],(*y)[i],(*z)[i]);
      }
    }
  }
}

//----------------------------------------------------------------------

void InitialDefault::copy_values_ 
(
 const FieldDescr * field_descr,
 FieldBlock *       field_block,
 double *           value, 
 bool *             region,
 int                index_field,
 int nx, int ny, int nz
 ) throw()
{

  // Copy the floating-point values to the field where logical values are true

  void * field = field_block->field_unknowns(field_descr,index_field);

  // Determine allocated array size

  int mx=0,my=0,mz=0;

  mx = nx;
  my = ny;
  mz = nz;
  if (field_block->ghosts_allocated()){
    int gx,gy,gz;
    field_descr->ghosts(index_field,&gx,&gy,&gz);
    mx += 2*gx;
    my += 2*gy;
    mz += 2*gz;
  }

  // Copy evaluated values to field values

  precision_enum precision = field_descr->precision(index_field);
  switch (precision) {
  case precision_single:
    for (int iz = 0; iz<nz; iz++) {
      for (int iy = 0; iy<ny; iy++) {
	for (int ix = 0; ix<nx; ix++) {
	  int in=ix + nx*(iy + ny*iz);
	  int im=ix + mx*(iy + my*iz);
	  if (region[in]) ((float *) field)[im] = (float) value[in];
	}
      }
    }
    break;
  case precision_double:
    for (int iz = 0; iz<nz; iz++) {
      for (int iy = 0; iy<ny; iy++) {
	for (int ix = 0; ix<nx; ix++) {
	  int in=ix + nx*(iy + ny*iz);
	  int im=ix + mx*(iy + my*iz);
	  if (region[in]) ((double *) field)[im] = (double) value[in];
	}
      }
    }
    break;
  case precision_quadruple:
    for (int iz = 0; iz<nz; iz++) {
      for (int iy = 0; iy<ny; iy++) {
	for (int ix = 0; ix<nx; ix++) {
	  int in=ix + nx*(iy + ny*iz);
	  int im=ix + mx*(iy + my*iz);
	  if (region[in]) ((long double *) field)[im] = (long double) value[in];
	}
      }
    }
    break;
  default:
    break;
  }
}

//----------------------------------------------------------------------

void InitialDefault::evaluate_float_ 
(FieldBlock * field_block, int index_list, std::string field_name, 
 int n, double * value, double * deflt,
 double * x, double * y, double * z, double * t) throw ()
{

  // parameter: Initial : <field> : value

  parameter_enum value_type = 
    parameters_->list_type(index_list,"value");

  if (value_type != parameter_float_expr &&
      value_type != parameter_float) {
	  	      
    ERROR1("InitialDefault::evaluate_float_", 
	   "Odd-index elements of %s must be floating-point expressions",
	   field_name.c_str());
  }

  // Evaluate the floating-point expression

  if (value_type == parameter_float) {
    double v = parameters_->list_value_float(index_list,"value",0.0);
    for (int i=0; i<n; i++) value[i] = v;
  } else {
    parameters_->list_evaluate_float
      (index_list,"value",n,value,deflt,x,y,z,t);
  }
}

//----------------------------------------------------------------------

void InitialDefault::evaluate_logical_ 
(FieldBlock * field_block, int index_list, std::string field_name, 
 int n, bool * value, bool * deflt,
 double * x, double * y, double * z, double * t) throw ()
{

  // parameter: Initial : <field> : value

  parameter_enum value_type = 
    parameters_->list_type(index_list,"value");

  if (value_type != parameter_logical_expr &&
      value_type != parameter_logical) {

    ERROR1("InitialDefault::evaluate_logical",
	   "Even-index elements of %s must be logical expressions",
	   field_name.c_str());
  }

  // Evaluate the logical expression

  if (value_type == parameter_logical) {
    double v = parameters_->list_value_logical(index_list,"value",0.0);
    for (int i=0; i<n; i++) value[i] = v;
  } else {
    parameters_->list_evaluate_logical
      (index_list,"value",n,value,deflt,x,y,z,t);
  }

}
