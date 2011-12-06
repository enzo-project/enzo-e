// See LICENSE_CELLO file for license and copyright information

/// @file     method_InitialDefault.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-02-16
/// @bug      Image mask blank for krummhorn: libpng version? (1.2.44 vs 1.2.46); note out.png is created
/// @todo     Parse initialization parameters once per Simulation rather than once per Block
/// @brief    Implementation of the InitialDefault class

#include "cello.hpp"

#include "method.hpp"

//----------------------------------------------------------------------

InitialDefault::InitialDefault(Parameters * parameters) throw ()
  : Initial (),
    parameters_(parameters)
{
}

//----------------------------------------------------------------------

void InitialDefault::enforce
(
 const Hierarchy  * hierarchy,
 const FieldDescr * field_descr,
 Block            * block
 ) throw()
{
  // Initialize Fields according to parameters

  //--------------------------------------------------
  parameters_->group_set(0,"Initial");
  //--------------------------------------------------

  FieldBlock *       field_block = block->field_block();

  double *value=0, *vdeflt=0, *x=0, *y=0, *z=0, *t=0;
  bool * mask=0, *rdeflt=0;
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

      ASSERT1("InitialDefault::enforce",
	     "Length of list parameter Initial:%s:value must be odd",
	     field_name.c_str(),
	     (list_length % 2) == 1);

      // Allocate arrays if needed
      if (value == NULL) {
	allocate_xyzt_(block,index_field,
		       &nx,&ny,&nz,
		       &value, &vdeflt,
		       &mask,&rdeflt,
		       &x,&y,&z,&t);
      }

      // Evaluate last non-conditional equation in list

      n = nx*ny*nz;

      evaluate_float_ (field_block, list_length-1, field_name, 
			n, value,vdeflt,x,y,z,t);

      for (int i=0; i<n; i++) mask[i] = true;

      copy_values_ (field_descr,field_block,value, mask,index_field,nx,ny,nz);

      // Evaluate conditional equations in list

      for (int index_list=0; index_list < list_length-1; index_list+=2) {

	evaluate_float_ (field_block, index_list, field_name, 
			  n, value,vdeflt,x,y,z,t);

	evaluate_logical_ (hierarchy,block,field_block, index_list + 1, field_name, 
			  n, mask,rdeflt,x,y,z,t);

	copy_values_ (field_descr,field_block,value, mask,index_field,nx,ny,nz);

      }

    } else if (parameter_type == parameter_unknown) {
      WARNING1("InitialDefault::enforce",  
	       "Uninitialized field %s",
	       field_name.c_str());
      
    } else {
      ERROR2("InitialDefault::enforce",
	     "Illegal parameter type %s when initializing field %s",
	     parameter_type_name[parameter_type],field_name.c_str());
    }
  }
  // Deallocate arrays if needed
  if (value != NULL) {
    delete [] value;
    delete [] vdeflt;
    delete [] mask;
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
 bool   ** mask,bool   ** rdeflt,
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
  (*mask)   = new bool [n];
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
	(*mask)[i] = false;
	(*rdeflt)[i]  = false;
	(*x)[i] = xm + ix*hx;
	(*y)[i] = ym + iy*hy;
	(*z)[i] = zm + iz*hz;
	(*t)[i] = 0.0;
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
 bool *             mask,
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
	  if (mask[in]) ((float *) field)[im] = (float) value[in];
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
	  if (mask[in]) ((double *) field)[im] = (double) value[in];
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
	  if (mask[in]) ((long double *) field)[im] = (long double) value[in];
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
(const Hierarchy * hierarchy,
 const Block * block,
 FieldBlock * field_block, int index_list, std::string field_name, 
 int n, bool * mask, bool * deflt,
 double * x, double * y, double * z, double * t) throw ()
{

  // parameter: Initial : <field> : value

  parameter_enum value_type = 
    parameters_->list_type(index_list,"value");

  bool v;
  const char * file;
  int nxb,nyb,nzb;
  int nx,ny;

  switch (value_type) {
  case parameter_logical:
    v = parameters_->list_value_logical(index_list,"value",false);
    for (int i=0; i<n; i++) mask[i] = v;
    break;
  case parameter_logical_expr:
    parameters_->list_evaluate_logical
      (index_list,"value",n,mask,deflt,x,y,z,t);
    break;
  case parameter_string:
    field_block->size(&nxb,&nyb,&nzb);
    ASSERT1("InitialDefault::evaluate_logical",
	   "mask file %s requires problem to be 2D",
	    field_name.c_str(),
	    nzb == 1);
    file = parameters_->list_value_string(index_list,"value","default");
    create_png_mask_(mask,hierarchy,block,file,nxb,nyb,&nx,&ny);
    break;
  default:
    ERROR3("InitialDefault::evaluate_logical",
	   "Even-index element %d of %s is of illegal type %d",
	   index_list,field_name.c_str(),value_type);
    break;
  }
}


void InitialDefault::create_png_mask_
(
 bool            * mask,
 const Hierarchy * hierarchy,
 const Block     * block,
 const char      * pngfile,
 int               nxb,
 int               nyb,
 int             * nx,
 int             * ny
 )
{
  pngwriter png;

  // Open the PNG file

  png.readfromfile(pngfile);

  // Get the PNG file size

  (*nx) = png.getwidth();
  (*ny) = png.getheight();

  // Clear the block mask

  int nb = nxb*nyb;
  for (int i=0; i<nb; i++) mask[i] = false;

  // Get the hierarchy's lower and upper extents

  double lower_h[2];
  hierarchy->lower(&lower_h[0],&lower_h[1],&lower_h[2]);
  double upper_h[2];
  hierarchy->upper(&upper_h[0],&upper_h[1],&upper_h[2]);

  // Get the hierarchy's size

  double size_h[2];
  size_h[0] = upper_h[0]-lower_h[0];
  size_h[1] = upper_h[1]-lower_h[1];

  // Get the block's lower and upper extents

  double lower_b[2];
  block->lower(&lower_b[0],&lower_b[1],&lower_b[2]);

  double upper_b[2];
  block->upper(&upper_b[0],&upper_b[1],&upper_b[2]);

  // get the offset between the block and the hierarchy

  double offset_b[2];
  offset_b[0] = lower_b[0]-lower_h[0];
  offset_b[1] = lower_b[1]-lower_h[1];

  // get the block's cell width

  double hb[2];
  block->field_block()->cell_width(block,&hb[0],&hb[1]);

  // Compute the block mask from the image pixel values

  for (int iy_b=0; iy_b<nyb; iy_b++) {
    int iy_h = int((*ny)*(iy_b*hb[1]+offset_b[1])/(size_h[1]));
    for (int ix_b=0; ix_b<nxb; ix_b++) {
      int ix_h = int((*nx)*(ix_b*hb[0]+offset_b[0])/(size_h[0]));

      int i_b = ix_b + nxb*iy_b;

      int r = png.read(ix_h+1,iy_h+1,1);
      int g = png.read(ix_h+1,iy_h+1,2);
      int b = png.read(ix_h+1,iy_h+1,3);

      mask[i_b] = (r+g+b > 0);

    }
  }

  png.close();
  return;
}
