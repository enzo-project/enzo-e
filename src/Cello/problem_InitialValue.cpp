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
 int cycle, double time) throw ()
  : Initial (cycle, time),
    parameters_(parameters)
{
  parameters_->group_set(0,"Initial");
  parameters_->group_set(1,"value");
}

//----------------------------------------------------------------------

void InitialValue::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  bool up = p.isUnpacking();

  // I'm not clear why the following does not create some kind of memory leak
  if (up) parameters_ = new Parameters;
  p | *parameters_;
}

//----------------------------------------------------------------------

void InitialValue::enforce_block ( Block * block,
				   const Hierarchy  * hierarchy ) throw()
{

  Initial::enforce_block(block,hierarchy);
  
  // Initialize Fields according to parameters

  ASSERT("InitialValue::enforce_block",
	 "Block does not exist",
	 block != NULL);
  
  //--------------------------------------------------
  parameters_->group_set(0,"Initial");
  parameters_->group_set(1,"value");
  //--------------------------------------------------

  Data *data = block->data();

  FieldData *       field_data = data->field_data();

  double *val_array = NULL;
  double *xc=NULL, *yc=NULL, *zc=NULL, *xf=NULL, *yf=NULL, *zf=NULL; 
  double t = block->time();

  Field field = data->field();
  int nx, ny, nz; // number of cells per axis in the active zone
  field.size(&nx,&ny,&nz);
  int gx, gy, gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  FieldDescr * field_descr = cello::field_descr();
  for (int index_field = 0;
       index_field < field_descr->field_count();
       index_field++) {
    
    std::string field_name = field_descr->field_name(index_field);

    // If Initial:<field_name>:value is a list, try parsing it

    // parameter: Initial : <field> : value

    parameter_type parameter_type = parameters_->type(field_name);

    if (parameter_type == parameter_float) {

      field_data->clear(field_descr,parameters_->value_float(field_name,0.0), 
			 index_field);

    } else if (parameter_type == parameter_list) {

      Value value(parameters_, field_name);

      int cx,cy,cz;
      field.centering(index_field, &cx,&cy,&cz);

      if (val_array == NULL){
        xc = new double [nx+2*gx];
        yc = new double [ny+2*gy];
        zc = new double [nz+2*gz];
        data->field_cells (xc, yc, zc, gx, gy, gz);

	// set val_array as array of 0s big enough for corner centered fields
        val_array = new double[(nx+2*gx+1)*(ny+2*gy+1)*(nz+2*gy+1)]();
      }

      if (xf == NULL && (cx!=0 || cy!=0 ||cz != 0) ){
        xf = new double [nx+2*gx+1];
        yf = new double [ny+2*gy+1];
        zf = new double [nz+2*gz+1];
        data->field_cell_faces (xf, yf, zf, gx, gy, gz);
      }

      double *x = (cx == 0) ? xc : xf;
      double *y = (cy == 0) ? yc : yf;
      double *z = (cz == 0) ? zc : zf;

      int ndx=nx+2*gx+cx;
      int ndy=ny+2*gy+cy;
      int ndz=nz+2*gz+cz;

      // Following convention of earlier version: initializing values in a
      // temporary array of doubles. Then the values are copied into the field
      // and casted to the appropriate value.

      // The cast to double * in the following line is redundant
      value.evaluate((double *)val_array, t,
		     ndx,ndx,x, 
		     ndy,ndy,y,
		     ndz,ndz,z);

      copy_values_(field_data,val_array,index_field,ndx,ndy,ndz);

    } else if (parameter_type == parameter_unknown) {
      if (block->index().is_root()) {
	WARNING1("InitialValue::enforce_block",  
		 "Uninitialized field %s",
		 field_name.c_str());
      }
    } else {
      ERROR2("InitialValue::enforce_block",
	     "Illegal parameter type %s when initializing field %s",
	     parameter_type_name[parameter_type],field_name.c_str());
    }
  }
  // Deallocate arrays if needed
  if (val_array != NULL) {
    delete [] val_array;
    delete [] xc;
    delete [] yc;
    delete [] zc;
  }
  if (xf != NULL) {
    delete [] xf;
    delete [] yf;
    delete [] zf;
  }
}

//======================================================================

void InitialValue::copy_values_ 
(
 FieldData *       field_data,
 double *           value,
 int                index_field,
 int nx, int ny, int nz
 ) throw()
{

  FieldDescr * field_descr = cello::field_descr();

  // Copy the floating-point values to the field where mask values are true

  void * array = field_data->unknowns(field_descr,index_field);

  // Determine allocated array size

  int gx,gy,gz;

  field_descr->ghost_depth(index_field,&gx,&gy,&gz);

  int offset = gx + nx*(gy + ny*gz);

  // Copy evaluated values to field values

  //  @@@@ BUG: PNG image input
  //  @@@@ IC's SMALLER BUT CENTERED IF OFFSET SET TO 0
  //  @@@@ IC's SHIFTED TO LOWER LEFT
  //  mask[0][0],[1][1],[2][2] == 1

  // The above comment is from an earlier version

  precision_type precision = field_descr->precision(index_field);
  switch (precision) {
  case precision_single:
    copy_precision_((float *)array,offset,value,nx,ny,nz);
    break;
  case precision_double:
    copy_precision_((double *)array,offset,value,nx,ny,nz);
    break;
  case precision_quadruple:
    copy_precision_((long double *)array,offset,value,nx,ny,nz);
    break;
  default:
    break;
  }
}

//----------------------------------------------------------------------

template<class T>
void InitialValue::copy_precision_
(T * array,
 int offset,
 double * value,
 int nx, int ny, int nz)
{
  for (int iz = 0; iz<nz; iz++) {
    for (int iy = 0; iy<ny; iy++) {
      for (int ix = 0; ix<nx; ix++) {
	int i = ix + nx*(iy + ny*iz);
	(array - offset)[i] = (T) value[i];
      }
    }
  }
}
