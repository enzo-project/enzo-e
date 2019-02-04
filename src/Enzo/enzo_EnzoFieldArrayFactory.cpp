#include "cello.hpp"
#include "enzo.hpp"

// There is redundancy, this can be condensed!

EFlt3DArray EnzoFieldArrayFactory::from_name(std::string field_name)
{
  Field field = block_->data()->field();
  const int id = field.field_id(field_name);
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);
  return EFlt3DArray((enzo_float *) field.values(field_name), mz, my, mx);
}

EFlt3DArray EnzoFieldArrayFactory::from_grouping(Grouping &grouping,
						 std::string group_name,
						 int index)
{
  Field field = block_->data()->field();
  int size = grouping.size(group_name);

  ASSERT1("EnzoFieldArrayFactory",
	  "\"%s\" is not the name of a real group\n",
	  group_name.c_str(), (size != 0));
  ASSERT3("EnzoFieldArrayFactory",
	  "index=%d is larger than %d, the size of group \"%s\"\n",
	  index, size, group_name.c_str(),(size>index));

  return from_name(grouping.item(group_name,index));
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

EFlt3DArray EnzoFieldArrayFactory::load_interior_bfieldi_field(
      Grouping &grouping, int dim)
{
  Field field = block_->data()->field();
  int size = grouping.size("bfield");
  ASSERT("EnzoFieldArrayFactory",
	 "grouping must contain a bfield group.",
	 size>0);
  ASSERT2("EnzoFieldArrayFactory",
	  "dim=%d is larger than %d, the size of group \"bfield\"\n",
	  dim, size,(size>dim));

  int dix = (dim == 0) ? 1 : 0;
  int diy = (dim == 1) ? 1 : 0;
  int diz = (dix == diy) ? 1 : 0;

  std::string field_name = grouping.item("bfield",dim);
  enzo_float *data =(enzo_float *)field.values(field_name);

  int mx,my,mz;
  cell_centered_dim_(field, field.field_id(field_name), &mx, &my, &mz);

  EFlt3DArray array;

  if (field.is_permanent(field.field_id(field_name))){
    EFlt3DArray temp_array(data,mz+diz,my+diy,mx+dix);
    array = temp_array.subarray(diz, mz, diy, my, dix, mx);
  } else {
    array = EFlt3DArray(data, mz-diz, my-diy, mx-dix);
  }
  return array;
}
