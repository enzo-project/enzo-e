#include "cello.hpp"
#include "enzo.hpp"

EFlt3DArray EnzoFieldArrayFactory::from_name(std::string field_name)
{
  Field field = block_->data()->field();
  const int id = field.field_id(field_name);
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);
  return EFlt3DArray((enzo_float *) field.values(field_name), mz, my, mx);
}

void check_grouping_details_(Grouping &grouping, std::string group_name,
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

EFlt3DArray EnzoFieldArrayFactory::from_grouping(Grouping &grouping,
						 std::string group_name,
						 int index)
{
  check_grouping_details_(grouping, group_name, index);
  return from_name(grouping.item(group_name,index));
}

EFlt3DArray EnzoFieldArrayFactory::reconstructed_field(Grouping &grouping,
						       std::string group_name,
						       int index, int dim)
{
  check_grouping_details_(grouping, group_name, index);
  ASSERT("EnzoFieldArrayFactory",
	 "reconstructed fields can only be face-centered along dim = 0, 1 or 2",
	 dim>-1 && dim <3);

  Field field = block_->data()->field();
  std::string field_name = grouping.item(group_name,index);
  const int id = field.field_id(field_name);
  int mx, my, mz;
  field.dimensions(id,&mx,&my,&mz);

  if (dim == 0){
    mx++;
  } else if (dim == 1){
    my++;
  } else {
    mz++;
  }
  return EFlt3DArray((enzo_float *) field.values(field_name), mz, my, mx);
}

EFlt3DArray EnzoFieldArrayFactory::interior_bfieldi(Grouping &grouping, int dim)
{
  EFlt3DArray t = from_grouping(grouping,"bfield",dim);
  if (dim == 0){
    return t.subarray(0, t.dim_size(2), 0, t.dim_size(1), 1, t.dim_size(0) - 1);
  } else if (dim == 1){
    return t.subarray(0, t.dim_size(2), 1, t.dim_size(1) - 1, 0, t.dim_size(0));
  } else {
    return t.subarray(1, t.dim_size(2) - 1, 1, t.dim_size(1), 0, t.dim_size(0));
  }
}

// Helper function that returns the dimensions of the cell-centered grid
// This doesn't take dimensions from EnzoBlock because the functions may be
// called while setting initial values (when EnzoBlock has not yet been set up)
void cell_centered_dim_(Field &field, int id, int *mx, int *my, int *mz)
{
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(id,&gx,&gy,&gz);

  *mx = nx + 2*gx;
  *my = ny + 2*gx;
  *mz = nz + 2*gx;
}

EFlt3DArray  EnzoFieldArrayFactory::load_temp_interface_grouping_field(
	Grouping &grouping, std::string group_name, int index,
	bool cell_centered_x, bool cell_centered_y, bool cell_centered_z)
{
  Field field = block_->data()->field();
  int size = grouping.size(group_name);
  ASSERT1("EnzoFieldArrayFactory",
	  "\"%s\" is not the name of a real group\n",
	  group_name.c_str(), (size != 0));
  ASSERT3("EnzoFieldArrayFactory",
	  "index=%d is larger than %d, the size of group \"%s\"\n",
	  index, size, group_name.c_str(),(size>index));

  std::string field_name = grouping.item(group_name,index);

  int mx, my, mz;
  // This is a very hacky short term solution to be used while making temporary
  // fields face-centered, again.
  cell_centered_dim_(field, field.field_id("density"), &mx, &my, &mz);
  if (!cell_centered_x) {
    mx--;
  }
  if (!cell_centered_y) {
    my--;
  }
  if (!cell_centered_z) {
    mz--;
  }
  return EFlt3DArray((enzo_float *)field.values(field_name),mz,my,mx);
}
