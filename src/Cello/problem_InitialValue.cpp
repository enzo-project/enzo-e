// See LICENSE_CELLO file for license and copyright information

/// @file     problem_InitialValue.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2011-02-16
/// @brief    Implementation of the InitialValue class

#include "cello.hpp"
#include <errno.h>
#include "problem.hpp"

//----------------------------------------------------------------------

InitialValue::InitialValue
(Parameters * parameters,
 int cycle, double time) throw()
    : Initial(cycle, time),
      ic_pairs_()
{
  parameters->group_set(0,"Initial");
  parameters->group_set(1,"value");

  std::vector<std::string> field_names = parameters->leaf_parameter_names();
  for (auto&& field_name : field_names){
    ic_pairs_.push_back
      ( std::make_pair(field_name, Value(parameters, field_name)) );
  }

}

//----------------------------------------------------------------------

void InitialValue::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change
  TRACEPUP;
  Initial::pup(p);
  p | ic_pairs_;
}

//----------------------------------------------------------------------

void InitialValue::enforce_block ( Block * block,
				   const Hierarchy  * hierarchy ) throw()
{

  Initial::enforce_block(block,hierarchy);

  ASSERT("InitialValue::enforce_block",
	 "Block does not exist",
	 block != NULL);

  Data *data = block->data();

  FieldData *       field_data = data->field_data();

  double *val_array = nullptr;
  double *xc=nullptr, *yc=nullptr, *zc=nullptr;
  double *xf=nullptr, *yf=nullptr, *zf=nullptr; 
  double t = block->time();

  Field field = data->field();
  int nx, ny, nz; // number of cells per axis in the active zone
  field.size(&nx,&ny,&nz);
  int gx, gy, gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  std::vector<bool> accessed_ic_pairs(ic_pairs_.size(), false);

  FieldDescr * field_descr = cello::field_descr();
  for (int index_field = 0;
       index_field < field_descr->field_count();
       index_field++) {
    
    std::string field_name = field_descr->field_name(index_field);

    // Search for the index of ic_pairs_ that holds the Value object that is to
    // be used to initialize the current field.
    std::size_t ic_value_index = ic_pairs_.size();
    for (std::size_t i = 0; i < ic_pairs_.size(); i++){
      if (ic_pairs_[i].first == field_name){
        ic_value_index = i;
        break;
      }
    }

    // Now, try to evaluate the initial conditions:
    if ((ic_value_index == ic_pairs_.size()) && (block->index().is_root())){
      // No initialization data was specified for the current field
      WARNING1("InitialValue::enforce_block", "Uninitialized field %s",
               field_name.c_str());

    } else if (ic_value_index < ic_pairs_.size()){
      // Initialization data was specified for the current field
      accessed_ic_pairs[ic_value_index] = true;
      const Value& cur_value = ic_pairs_[ic_value_index].second;

      if (cur_value.wraps_single_float_param()) {
	field_data->clear(field_descr, cur_value.evaluate(t, 0.0, 0.0, 0.0),
			  index_field);
      } else {
	int cx,cy,cz;
	field.centering(index_field, &cx,&cy,&cz);

	if (val_array == nullptr){
	  xc = new double [nx+2*gx];
	  yc = new double [ny+2*gy];
	  zc = new double [nz+2*gz];
	  data->field_cells (xc, yc, zc, gx, gy, gz);

	  // set val_array as array of 0s big enough for corner centered fields
	  val_array = new double[(nx+2*gx+1)*(ny+2*gy+1)*(nz+2*gz+1)]();
	}

	if (xf == nullptr && (cx!=0 || cy!=0 ||cz != 0) ){
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
	// temporary array of doubles. Then the values are copied into the
	// field and casted to the appropriate value.

	// The cast to double * in the following line is redundant
	cur_value.evaluate((double *)val_array, t,
                           ndx,ndx,x,  ndy,ndy,y,  ndz,ndz,z);

	copy_values_(field_data,val_array,index_field,ndx,ndy,ndz);
      } 
    }
  }
  // Deallocate arrays if needed
  if (val_array != nullptr) {
    delete [] val_array;
    delete [] xc;
    delete [] yc;
    delete [] zc;
  }
  if (xf != nullptr) {
    delete [] xf;
    delete [] yf;
    delete [] zf;
  }

  // Finaly, confirm initialization data isn't specified for missing fields.
  for (std::size_t i = 0; i < ic_pairs_.size(); i++){
    if (!accessed_ic_pairs[i]){
      ERROR1("InitialValue::enforce_block",
             ("Initial conditions were specified for \"%s\", but that field "
              "doesn't exist"),
             ic_pairs_[i].first.c_str());
    }
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
