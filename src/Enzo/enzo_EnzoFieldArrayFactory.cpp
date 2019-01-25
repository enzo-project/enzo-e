#include "cello.hpp"
#include "enzo.hpp"

// There is a lot of redundancy, this can all be condensed!

void EnzoFieldArrayFactory::initialize_field_array(Block *block,
						   EnzoArray<enzo_float> &array,
						   std::string field_name)
{
  Field field = block->data()->field();
  const int id = field.field_id(field_name);

  // get the field dimensions
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);
  enzo_float *data = (enzo_float *) field.values(field_name);
  array.initialize_wrapper(data, mz, my, mx);
}

void EnzoFieldArrayFactory::load_grouping_field(Block *block,
						Grouping &grouping,
						std::string group_name,
						int index,
						EnzoArray<enzo_float> &array)
{
  Field field = block->data()->field();
  int size = grouping.size(group_name);

  ASSERT1("EnzoFieldArrayFactory",
	  "\"%s\" is not the name of a real group\n",
	  group_name.c_str(), (size != 0));
  ASSERT3("EnzoFieldArrayFactory",
	  "index=%d is larger than %d, the size of group \"%s\"\n",
	  index, size, group_name.c_str(),(size>index));
  std::string field_name = grouping.item(group_name,index);
  
  enzo_float *data = (enzo_float *) field.values(field_name);
  int id = field.field_id(field_name);
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);

  array.initialize_wrapper(data,mz,my,mx);
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

void EnzoFieldArrayFactory::load_temp_interface_grouping_field(
	Block *block, Grouping &grouping, std::string group_name, int index,
	EnzoArray<enzo_float> &array, bool cell_centered_x,
	bool cell_centered_y, bool cell_centered_z)
{
  int size = grouping.size(group_name);
  ASSERT1("EnzoFieldArrayFactory",
	  "\"%s\" is not the name of a real group\n",
	  group_name.c_str(), (size != 0));
  ASSERT3("EnzoFieldArrayFactory",
	  "index=%d is larger than %d, the size of group \"%s\"\n",
	  index, size, group_name.c_str(),(size>index));

  Field field = block->data()->field();
  std::string field_name = grouping.item(group_name,index);
  enzo_float *data =(enzo_float *)field.values(field_name);

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
  array.initialize_wrapper(data,mz,my,mx);
}

void EnzoFieldArrayFactory::load_interior_bfieldi_field(
      Block *block, Grouping &grouping, int dim,
      EnzoArray<enzo_float> &array)
{
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

  Field field = block->data()->field();
  std::string field_name = grouping.item("bfield",dim);
  enzo_float *data =(enzo_float *)field.values(field_name);

  int mx,my,mz;
  cell_centered_dim_(field, field.field_id(field_name), &mx, &my, &mz);
  
  if (field.is_permanent(field.field_id(field_name))){
    EnzoArray<enzo_float> temp_array;
    temp_array.initialize_wrapper(data,mz+diz,my+diy,mx+dix);
    array.initialize_subarray(temp_array, diz, mz, diy, my, dix, mx);
  } else {
    array.initialize_wrapper(data, mz-diz, my-diy, mx-dix);
  }
}
