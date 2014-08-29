// See LICENSE_CELLO file for license and copyright information

/// @file     problem_InitialValue.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-02-16
/// @brief    Implementation of the InitialValue class

#include "cello.hpp"
#include <errno.h>
#include "problem.hpp"

//----------------------------------------------------------------------

InitialValue::InitialValue
(Parameters * parameters,
 const FieldDescr * field_descr,
 int cycle, double time) throw ()
  : Initial (cycle, time),
    parameters_(parameters),
    field_descr_(field_descr),
    num_fields_(field_descr->field_count()),
    num_masks_(NULL),
    mask_(0),
    nx_(0),
    ny_(0)
{
  parameters_->group_set(0,"Initial");

  mask_      = new bool ** [num_fields_];
  nx_        = new int * [num_fields_];
  ny_        = new int * [num_fields_];
  num_masks_ = new int [num_fields_];
  for (int i=0; i<num_fields_; i++) {
    mask_[i] = 0;
    nx_[i] = 0;
    ny_[i] = 0;
    num_masks_[i] = 0;
  }
  
  for (int index_field = 0; index_field < num_fields_; index_field++) {

    std::string field_name = field_descr->field_name(index_field);

    parameters_->group_set(1,field_name);

    if (parameters_->type("value") == parameter_list) {
      int num_values = parameters_->list_length("value");
      if (num_values > 1) {
	num_masks_[index_field] =  num_values / 2;
	mask_[index_field] = new bool * [num_masks_[index_field]];
	nx_[index_field]   = new int [num_masks_[index_field]];
	ny_[index_field]   = new int [num_masks_[index_field]];
	// loop through value masks
	for (int index_mask=0; 
	     index_mask < num_masks_[index_field];
	     index_mask++) {
	  int index_value = index_mask*2+1;
	  if (parameters_->list_type(index_value,"value") == parameter_string) {
	    std::string file 
	      = parameters_->list_value_string(index_value,"value","default");
	    create_mask_png_ (&mask_[index_field][index_mask],
			      &nx_[index_field][index_mask],
			      &ny_[index_field][index_mask],
			      file);
			   
	  
	  }
	}
      }
    }
  }
}

//----------------------------------------------------------------------

InitialValue::~InitialValue() throw()
{
  for (int index_field = 0; index_field<num_fields_; index_field++) {
    for (int index_mask = 0; index_mask<num_masks_[index_field]; index_mask++) {
      if (mask_[index_field][index_mask]) {
	delete [] mask_[index_field][index_mask];
	mask_[index_field][index_mask] = 0;
      }
    }
    delete [] mask_[index_field]; mask_[index_field] = 0;
    delete [] nx_[index_field];   nx_[index_field]   = 0;
    delete [] ny_[index_field];   ny_[index_field]   = 0;
  }
  delete [] mask_; mask_ = 0;
  delete [] nx_;   nx_   = 0;
  delete [] ny_;   ny_   = 0;
}
//----------------------------------------------------------------------

void InitialValue::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  bool up = p.isUnpacking();

  if (up) parameters_ = new Parameters;
  p | *parameters_;
  if (up) field_descr_ = new FieldDescr;
  p | *((FieldDescr *) field_descr_);
  p | num_fields_;
  WARNING("InitialValue::pup","mask_[][] not pupped");
  if (up) num_masks_ = new int[num_fields_];
  PUParray(p,num_masks_,num_fields_);
}

//----------------------------------------------------------------------

void InitialValue::enforce_block
(
 CommBlock        * comm_block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()
{
  // Initialize Fields according to parameters

  ASSERT("InitialValue::enforce_block",
	 "CommBlock does not exist",
	 comm_block != NULL);
  
  //--------------------------------------------------
  parameters_->group_set(0,"Initial");
  //--------------------------------------------------

  FieldBlock *       field_block = comm_block->block()->field_block();

  double *value=0, *vdeflt=0, *x=0, *y=0, *z=0, t;
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

    parameter_type parameter_type = parameters_->type("value");

    if (parameter_type == parameter_float) {

      field_block->clear(parameters_->value_float("value",0.0), 
			 index_field);

    } else if (parameter_type == parameter_list) {

      // Check parameter length

      int list_length = parameters_->list_length("value");

      ASSERT1("InitialValue::enforce_block",
	     "Length of list parameter Initial:%s:value must be odd",
	     field_name.c_str(),
	     (list_length % 2) == 1);

      // Allocate arrays if needed
      if (value == NULL) {
	allocate_xyzt_(comm_block,index_field,
		       field_block,
		       field_descr,
		       &nx,&ny,&nz,
		       &value, &vdeflt,
		       &mask,  &rdeflt,
		       &x,&y,&z,&t);
      }

      // Evaluate last non-conditional equation in list

      n = nx*ny*nz;

      evaluate_float_ (field_block, list_length-1, field_name, 
			n, value,vdeflt,x,y,z,t);

      copy_values_ (field_descr,field_block,value, NULL,index_field,nx,ny,nz);

      // Evaluate conditional equations in list

      for (int index_value=0; index_value < list_length-1; index_value+=2) {

	evaluate_float_ (field_block, index_value, field_name, 
			  n, value,vdeflt,x,y,z,t);

	evaluate_mask_ 
	  (hierarchy,comm_block,field_block, index_field, index_value+1,
	   field_name,   field_descr, n, mask,rdeflt,x,y,z,t);

	copy_values_ (field_descr,field_block,value, mask,index_field,nx,ny,nz);

      }

    } else if (parameter_type == parameter_unknown) {
      static bool displayed[MAX_FIELDS] = {false};
      if (! displayed[index_field]) {
      WARNING1("InitialValue::enforce_block",  
	       "Uninitialized field %s",
	       field_name.c_str());
      displayed[index_field] = true;
      }
      
    } else {
      ERROR2("InitialValue::enforce_block",
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
  }
}

//======================================================================

void InitialValue::allocate_xyzt_
(
 CommBlock * comm_block,
 int index_field,
 const FieldBlock * field_block,
 const FieldDescr * field_descr,
 int * nx, int * ny, int * nz,
 double ** value, double ** vdeflt,
 bool   ** mask,bool   ** rdeflt,
 double ** x, double ** y, double ** z, double * t
 ) throw()
{

  // Get field size

  field_block->size(nx,ny,nz);

  int gx,gy,gz;
  field_descr->ghosts(index_field,&gx,&gy,&gz);
  (*nx) += 2*gx;
  (*ny) += 2*gy;
  (*nz) += 2*gz;

  int n = (*nx)*(*ny)*(*nz);

  // Allocate double arrays

  (*value)  = new double [n];
  (*vdeflt) = new double [n];
  (*mask)   = new bool [n];
  (*rdeflt) = new bool [n];
  (*x)      = new double [n];
  (*y)      = new double [n];
  (*z)      = new double [n];

  double xm, xp, ym, yp, zm, zp;

  comm_block->block()->lower(&xm,&ym,&zm);
  comm_block->block()->upper(&xp,&yp,&zp);

  double hx,hy,hz;
  field_block->cell_width(xm,xp,&hx,
			  ym,yp,&hy,
			  zm,zp,&hz);

  *t = comm_block->time();

  // Initialize arrays
  for (int iz=0; iz<(*nz); iz++) {
    for (int iy=0; iy<(*ny); iy++) {
      for (int ix=0; ix<(*nx); ix++) {
	int i=ix + (*nx)*(iy + (*ny)*iz);
	(*value)[i]  = 0.0;
	(*vdeflt)[i] = 0.0;
	(*mask)[i]   = false;
	(*rdeflt)[i] = false;
	(*x)[i]      = xm + (ix-gx+0.5)*hx;
	(*y)[i]      = ym + (iy-gy+0.5)*hy;
	(*z)[i]      = zm + (iz-gz+0.5)*hz;
      }
    }
  }
}

//----------------------------------------------------------------------

void InitialValue::copy_values_ 
(
 const FieldDescr * field_descr,
 FieldBlock *       field_block,
 double *           value, 
 bool *             mask,
 int                index_field,
 int nx, int ny, int nz
 ) throw()
{

  // Copy the floating-point values to the field where mask values are true

  void * field = field_block->unknowns(index_field);

  // Determine allocated array size

  int gx,gy,gz;

  field_descr->ghosts(index_field,&gx,&gy,&gz);

  int offset = gx + nx*(gy + ny*gz);

  // Copy evaluated values to field values

  //  @@@@ BUG: PNG image input
  //  @@@@ IC's SMALLER BUT CENTERED IF OFFSET SET TO 0
  //  @@@@ IC's SHIFTED TO LOWER LEFT
  //  mask[0][0],[1][1],[2][2] == 1

  precision_type precision = field_descr->precision(index_field);
  switch (precision) {
  case precision_single:
    if (mask) copy_precision_((float *)field,mask,offset,value,nx,ny,nz);
    else      copy_precision_((float *)field,     offset,value,nx,ny,nz);
    break;
  case precision_double:
    if (mask) copy_precision_((double *)field,mask,offset,value,nx,ny,nz);
    else      copy_precision_((double *)field,     offset,value,nx,ny,nz);
    break;
  case precision_quadruple:
    if (mask) copy_precision_((long double *)field,mask,offset,value,nx,ny,nz);
    else      copy_precision_((long double *)field,     offset,value,nx,ny,nz);
    break;
  default:
    break;
  }
}

//----------------------------------------------------------------------

template<class T>
void InitialValue::copy_precision_
(T * field,
 bool * mask,
 int offset,
 double * value,
 int nx, int ny, int nz)
{
  for (int iz = 0; iz<nz; iz++) {
    for (int iy = 0; iy<ny; iy++) {
      for (int ix = 0; ix<nx; ix++) {
	int i = ix + nx*(iy + ny*iz);
	if (mask[i]) 
	  (field - offset)[i] = (T) value[i];
      }
    }
  }
}

//----------------------------------------------------------------------

//----------------------------------------------------------------------

template<class T>
void InitialValue::copy_precision_
(T * field,
 int offset,
 double * value,
 int nx, int ny, int nz)
{
  for (int iz = 0; iz<nz; iz++) {
    for (int iy = 0; iy<ny; iy++) {
      for (int ix = 0; ix<nx; ix++) {
	int i = ix + nx*(iy + ny*iz);
	(field - offset)[i] = (T) value[i];
      }
    }
  }
}

//----------------------------------------------------------------------

void InitialValue::evaluate_float_ 
(FieldBlock * field_block, int index_value, std::string field_name, 
 int n, double * value, double * deflt,
 double * x, double * y, double * z, double t) throw ()
{

  // parameter: Initial : <field> : value

  parameter_type value_type = 
    parameters_->list_type(index_value,"value");

  if (value_type != parameter_float_expr &&
      value_type != parameter_float) {
	  	      
    ERROR1("InitialValue::evaluate_float_", 
	   "Odd-index elements of %s must be floating-point expressions",
	   field_name.c_str());
  }

  // Evaluate the floating-point expression

  if (value_type == parameter_float) {
    double v = parameters_->list_value_float(index_value,"value",0.0);
    for (int i=0; i<n; i++) value[i] = v;
  } else {
    parameters_->list_evaluate_float
      (index_value,"value",n,value,deflt,x,y,z,t);
  }
}

//----------------------------------------------------------------------

void InitialValue::evaluate_mask_ 
(const Hierarchy * hierarchy,
 const CommBlock * comm_block,
 FieldBlock * field_block, 
 int index_field,  int index_value,
 std::string field_name, 
 const FieldDescr * field_descr,
 int n, bool * mask, bool * deflt,
 double * x, double * y, double * z, double t) throw ()
{

  // parameter: Initial : <field> : value

  int index_mask = index_value / 2;

  parameter_type value_type = 
    parameters_->list_type(index_value,"value");

  bool v;
  std::string file;
  int nxb,nyb,nzb;

  switch (value_type) {
  case parameter_logical:

    // logical constant

    v = parameters_->list_value_logical(index_value,"value",false);
    for (int i=0; i<n; i++) mask[i] = v;
    break;

  case parameter_logical_expr:

    // logical expression

    parameters_->list_evaluate_logical
      (index_value,"value",n,mask,deflt,x,y,z,t);
    break;

  case parameter_string:

    {
      field_block->size(&nxb,&nyb,&nzb);
      ASSERT1("InitialValue::evaluate_logical",
	      "mask file %s requires problem to be 2D",
	      field_name.c_str(),
	      nyb > 1 && nzb == 1);

      bool * mask_png = mask_[index_field][index_mask];
      int nx_png = nx_[index_field][index_mask];
      int ny_png = ny_[index_field][index_mask];
      evaluate_mask_png_
	(mask,nxb,nyb,
	 mask_png,nx_png, ny_png,
	 hierarchy,comm_block,field_descr);
    }
    break;

  default:
    ERROR3("InitialValue::evaluate_mask",
	   "Even-index element %d of %s is of illegal type %d",
	   index_value,field_name.c_str(),value_type);
    break;
  }
}

//----------------------------------------------------------------------

void InitialValue::evaluate_mask_png_
(
 bool            * mask_block, int nxb, int nyb,
 bool            * mask_png,   int nx,  int ny,
 const Hierarchy * hierarchy,
 const CommBlock * comm_block,
 const FieldDescr * field_descr
 )
{
  int gx,gy,gz;
  field_descr->ghosts(0,&gx,&gy,&gz);
  nxb += 2*gx;
  nyb += 2*gy;

  // Clear the block mask

  int nb = nxb*nyb;
  for (int i=0; i<nb; i++) mask_block[i] = false;

  // Get the hierarchy's lower and upper extents

  double lower_h[3];
  hierarchy->lower(&lower_h[0],&lower_h[1],&lower_h[2]);
  double upper_h[3];
  hierarchy->upper(&upper_h[0],&upper_h[1],&upper_h[2]);

  // Get the block's lower and upper extents

  double lower_b[3];
  comm_block->block()->lower(&lower_b[0],&lower_b[1],&lower_b[2]);

  double upper_b[3];
  comm_block->block()->upper(&upper_b[0],&upper_b[1],&upper_b[2]);

  // get the block's cell width

  double hb[3];
  comm_block->block()->field_block()->cell_width 
    (lower_b[0],upper_b[0],&hb[0],
     lower_b[1],upper_b[1],&hb[1],
     lower_b[2],upper_b[2],&hb[2]);

  // Get the hierarchy's size including ghosts

  double size_h[3];
  size_h[0] = upper_h[0]-lower_h[0] + 2*gx*hb[0];
  size_h[1] = upper_h[1]-lower_h[1] + 2*gy*hb[1];

  // get the offset between the block and the hierarchy

  double offset_b[3];
  offset_b[0] = lower_b[0]-lower_h[0];
  offset_b[1] = lower_b[1]-lower_h[1];

  for (int iy_b=0; iy_b<nyb; iy_b++) {
    int iy_h = int(ny*(iy_b*hb[1]+offset_b[1])/(size_h[1]));
    for (int ix_b=0; ix_b<nxb; ix_b++) {
      int ix_h = int(nx*(ix_b*hb[0]+offset_b[0])/(size_h[0]));
      
      int i_b = ix_b + nxb*(iy_b);
      int i_h = ix_h + nx*iy_h;

      mask_block[i_b] = mask_png[i_h];
    }
  }
  return;
}

//----------------------------------------------------------------------

void InitialValue::create_mask_png_
( bool ** mask,  int * nx, int * ny, std::string pngfile)
{
  pngwriter png;

  // Open the PNG file

  errno = 0;
  png.readfromfile(pngfile.c_str());
  ASSERT1 ("InitialValue::create_mash_png_()",
	   "Cannot read PNG file %s",
	   pngfile.c_str(),
	   errno == 0);

  // Get the PNG file size

  (*nx) = png.getwidth();
  (*ny) = png.getheight();

  // Allocate and clear the mask

  int n = (*nx)*(*ny);
  (*mask) = new bool [n];
  for (int i=0; i<n; i++) (*mask)[i] = false;

  for (int iy=0; iy<(*ny); iy++) {
    for (int ix=0; ix<(*nx); ix++) {

      int i = ix + (*nx)*iy;

      int r = png.read(ix+1,iy+1,1);
      int g = png.read(ix+1,iy+1,2);
      int b = png.read(ix+1,iy+1,3);

      (*mask)[i] = (r+g+b > 0);
    }
  }
  png.close();
  return;
}
