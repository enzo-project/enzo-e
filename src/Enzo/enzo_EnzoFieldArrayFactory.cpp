// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFieldArrayFactory.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implementation of EnzoFieldArrayFactory class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::from_name(std::string field_name)
{
  EFlt3DArray temp = full_field_from_name_(field_name);
  if (stale_depth_ != 0){
    return exclude_stale_cells_(temp);
  } else {
    return temp;
  }
}

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::from_grouping(Grouping &grouping,
						 std::string group_name,
						 int index)
{
  check_grouping_details_(grouping, group_name, index);
  return from_name(grouping.item(group_name,index));
}

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::assigned_center_from_grouping
(Grouping &grouping, std::string group_name, int index, int dim)
{
  check_grouping_details_(grouping, group_name, index);
  ASSERT("EnzoFieldArrayFactory",
	 "reconstructed fields can only be face-centered along dim = 0, 1 or 2",
	 dim>-1 && dim <3);

  Field field = block_->data()->field();
  std::string field_name = grouping.item(group_name,index);
  const int id = field.field_id(field_name);

  int c_arr[3]; // = {cx, cy, cz}
  field.centering(id, c_arr, c_arr+1, c_arr+2);
  int temp_sum = c_arr[0]+c_arr[1]+c_arr[2];
  if ((c_arr[0] > 0) || (c_arr[1] > 0) || (c_arr[2] > 0)){

    ERROR1("EnzoFieldArrayFactory::assigned_center_from_grouping",
	   "\"%s\" has space for face-centered values on the grid exterior",
	   field_name.c_str());

  } else if (temp_sum < -1){

    ERROR1("EnzoFieldArrayFactory::assigned_center_from_grouping",
	   "the \"%s\" field is face-centered along multiple dimensions",
	   field_name.c_str());

  } else if ((temp_sum == -1) && (c_arr[dim] != -1)){

    int actual_face = -1;
    for (int i=0; i<3; i++){ if (c_arr[i] == -1) { actual_face = i; } }
    char axis_names[3] = {'x','y','z'};

    ERROR3("EnzoFieldArrayFactory::assigned_center_from_grouping",
	   ("the underlying \"%s\" field must be cell-centered or face-"
	    "centered along the %c-axis; not face-centered along the %c-axis"),
	   field_name.c_str(), axis_names[dim], axis_names[actual_face]);
  }

  // separate the following if-statement from the above to avoid compiler
  // warning about the return type

  if (temp_sum == -1){

    // the field's underlying centering already matches the specified centering
    return from_name(field_name);

  } else {

    // the underlying field is cell-centered. Modify it's shape to match the
    // designated centering.

    int mx,my,mz;
    field.dimensions(id,&mx,&my,&mz);
    if (dim == 0){
      mx--;
    } else if (dim == 1){
      my--;
    } else {
      mz--;
    }

    EFlt3DArray temp((enzo_float *) field.values(field_name), mz, my, mx);
    if (stale_depth_ != 0){
      return exclude_stale_cells_(temp);
    } else {
      return temp;
    }
  }
}

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::interior_bfieldi(Grouping &grouping, int dim)
{
  check_grouping_details_(grouping, "bfield", dim);
  EFlt3DArray t = full_field_from_name_(grouping.item("bfield",dim));

  if (dim == 0){
    return t.subarray(CSlice(stale_depth_,     t.shape(0) - stale_depth_),
		      CSlice(stale_depth_,     t.shape(1) - stale_depth_),
		      CSlice(1 + stale_depth_, t.shape(2) - 1 - stale_depth_));
  } else if (dim == 1){
    return t.subarray(CSlice(stale_depth_,     t.shape(0) - stale_depth_),
		      CSlice(1 + stale_depth_, t.shape(1) - 1 - stale_depth_),
		      CSlice(stale_depth_,     t.shape(2) - stale_depth_));
  } else {
    return t.subarray(CSlice(1 + stale_depth_, t.shape(0) - 1 - stale_depth_),
		      CSlice(stale_depth_,     t.shape(1) - stale_depth_),
		      CSlice(stale_depth_,     t.shape(2) - stale_depth_));
  }
}

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::full_field_from_name_(std::string field_name)
{
  Field field = block_->data()->field();
  const int id = field.field_id(field_name);
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);
  EFlt3DArray temp((enzo_float *) field.values(field_name), mz, my, mx);
  return temp;
}

//----------------------------------------------------------------------

void EnzoFieldArrayFactory::check_grouping_details_(Grouping &grouping,
						    std::string group_name,
						    int index)
{
  int size = grouping.size(group_name);
  ASSERT1("EnzoFieldArrayFactory",
	  "\"%s\" is not the name of a real group\n",
	  group_name.c_str(), (size != 0));
  ASSERT3("EnzoFieldArrayFactory",
	  "index=%d is larger than %d, the size of group \"%s\"\n",
	  index, size, group_name.c_str(),(size>index));
}

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::exclude_stale_cells_(EFlt3DArray &arr)
{
  ASSERT("EnzoFieldArrayFactory","each dim of arr must exceed 2*stale_depth_.",
	 2*stale_depth_ < arr.shape(0) && 2*stale_depth_ < arr.shape(1) &&
	 2*stale_depth_ < arr.shape(2));

  return arr.subarray(CSlice(stale_depth_, -stale_depth_),
		      CSlice(stale_depth_, -stale_depth_),
		      CSlice(stale_depth_, -stale_depth_));
}
