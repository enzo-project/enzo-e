// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFieldArrayFactory.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implementation of EnzoFieldArrayFactory class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::from_name(std::string field_name,
					     int stale_depth)
{
  Field field = block_->data()->field();
  const int id = field.field_id(field_name);
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);
  EFlt3DArray temp((enzo_float *) field.values(field_name), mz, my, mx);
  if (stale_depth != 0){
    return exclude_stale_cells_(temp,stale_depth);
  } else {
    return temp;
  }
}

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::from_grouping(Grouping &grouping,
						 std::string group_name,
						 int index, int stale_depth)
{
  check_grouping_details_(grouping, group_name, index);
  return from_name(grouping.item(group_name,index),stale_depth);
}

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::reconstructed_field(Grouping &grouping,
						       std::string group_name,
						       int index, int dim,
						       int stale_depth)
{
  check_grouping_details_(grouping, group_name, index);
  ASSERT("EnzoFieldArrayFactory",
	 "reconstructed fields can only be face-centered along dim = 0, 1 or 2",
	 dim>-1 && dim <3);

  Field field = block_->data()->field();
  std::string field_name = grouping.item(group_name,index);
  const int id = field.field_id(field_name);
  int mx, my, mz;
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(id,&gx,&gy,&gz);

  mx = nx + 2*gx;
  my = ny + 2*gy;
  mz = nz + 2*gz;
  

  if (dim == 0){
    mx--;
  } else if (dim == 1){
    my--;
  } else {
    mz--;
  }

  EFlt3DArray temp((enzo_float *) field.values(field_name), mz, my, mx);
  if (stale_depth != 0){
    return exclude_stale_cells_(temp,stale_depth);
  } else {
    return temp;
  }
}

//----------------------------------------------------------------------

EFlt3DArray EnzoFieldArrayFactory::interior_bfieldi(Grouping &grouping, int dim,
						    int stale_depth)
{
  EFlt3DArray t = from_grouping(grouping,"bfield",dim);
  if (dim == 0){
    return t.subarray(ESlice(stale_depth,     t.shape(0) - stale_depth),
		      ESlice(stale_depth,     t.shape(1) - stale_depth),
		      ESlice(1 + stale_depth, t.shape(2) - 1 - stale_depth));
  } else if (dim == 1){
    return t.subarray(ESlice(stale_depth,     t.shape(0) - stale_depth),
		      ESlice(1 + stale_depth, t.shape(1) - 1 - stale_depth),
		      ESlice(stale_depth,     t.shape(2) - stale_depth));
  } else {
    return t.subarray(ESlice(1 + stale_depth, t.shape(0) - 1 - stale_depth),
		      ESlice(stale_depth,     t.shape(1) - stale_depth),
		      ESlice(stale_depth,     t.shape(2) - stale_depth));
  }
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

EFlt3DArray EnzoFieldArrayFactory::exclude_stale_cells_(EFlt3DArray &arr,
							int stale_depth)
{
  ASSERT("EnzoFieldArrayFactory", "stale_depth must be a positive int",
	 stale_depth > 0);
  ASSERT("EnzoFieldArrayFactory", "each axis of arr must exceed 2*stale_depth.",
	 2*stale_depth < arr.shape(0) && 2*stale_depth < arr.shape(1) &&
	 2*stale_depth < arr.shape(2));

  return arr.subarray(ESlice(stale_depth, -stale_depth),
		      ESlice(stale_depth, -stale_depth),
		      ESlice(stale_depth, -stale_depth));
}
