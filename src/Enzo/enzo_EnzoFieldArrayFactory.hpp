// This will consolidate the various creational helper functions
// used to construct EnzoArrays that wrap field values.
// For now, this is purely organizational. But in the future:
//  - refactor to take a block or field in the constructor.
//  - refactor once we allow temporary fields to be face-centered again.
//  - Now that they are grouped under EnzoFieldArrayFactory, they can be
//    renamed with shorter names (which are equally as descriptive

#ifndef ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP
#define ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP

// alias for EnzoArray<enzo_float,3>
typedef EnzoArray<enzo_float,3> EFlt3DArray;


class EnzoFieldArrayFactory
{
  /// @class    EnzoFieldArrayFactory
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates construction of EnzoArrays from fields

public:
  EnzoFieldArrayFactory(Block *block)
  {
    block_ = block;
  }

  ~EnzoFieldArrayFactory()
  { block_ = NULL;}

  // get a field as an array
  EFlt3DArray from_name(std::string field_name);

  // reading in Grouping fields
  EFlt3DArray from_grouping(Grouping &grouping, std::string group_name,
			    int index);

  // This is only used to load in reconstructed arrays (name is a relic from
  // a period when I incorrectly believed that temporary arrays could not be
  // face-centered) Setting cell_centered_{dim} to true indicates that values
  // are cell-centered along {dim}
  // Going to rename this and simplify the function signature. For
  // reconstructed data, we only need to indicate the face-centered direction
  EFlt3DArray load_temp_interface_grouping_field(Grouping &grouping,
						 std::string group_name,
						 int index,
						 bool cell_centered_x,
						 bool cell_centered_y,
						 bool cell_centered_z);
  
  // Initialize an EnzoArray representing a component of the interface B-fields.
  // This function accepts both the grouping of permanent and temporary fields,
  // and in each case yields an array that does not include values for the
  // exterior face of the grid. For example, dim=0, the dimension of the
  // resulting EnzoArray is (mz, my-1, mx)
  EFlt3DArray load_interior_bfieldi_field(Grouping &grouping, int dim);

private:
  Block* block_;
};
#endif /* ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP */
