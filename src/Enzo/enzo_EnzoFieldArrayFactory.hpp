// This will consolidate the various creational helper functions
// used to construct EnzoArrays that wrap field values.
// For now, this is purely organizational. But in the future:
//  - refactor to take a block or field in the constructor.
//  - refactor once we allow temporary fields to be face-centered again.
//  - Now that they are grouped under EnzoFieldArrayFactory, they can be
//    renamed with shorter names (which are equally as descriptive

#ifndef ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP
#define ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP

// alias for EnzoArray<enzo_float>
typedef EnzoArray<enzo_float> EnzoFltArray;


class EnzoFieldArrayFactory
{
  /// @class    EnzoFieldArrayFactory
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates construction of EnzoArrays from fields

public:
  EnzoFieldArrayFactory() throw()
  {}

  // get a field as an array
  void initialize_field_array(Block *block, EnzoArray<enzo_float> &array,
			      std::string field_name);
  void initialize_field_array(Block *block, EnzoArray<enzo_float> &array,
			      int field_id);

  // reading in Grouping fields
  void load_grouping_field(Block *block, Grouping &grouping,
			   std::string group_name, int index,
			   EnzoArray<enzo_float> &array);

  // Because temporary fields must be allocated as cell-centered, we do not
  // have FieldDescr internally track the centering of the fields. Consequently,
  // we must specify the centering of the EnzoArray. Setting cell_centered_{dim}
  // to true, it indicates that values are cell-centered along {dim}. Face
  // -centered temporary fields do not have values at the exterior of the mesh.
  void load_temp_interface_grouping_field(Block *block, Grouping &grouping,
					  std::string group_name, int index,
					  EnzoArray<enzo_float> &array,
					  bool cell_centered_x,
					  bool cell_centered_y,
					  bool cell_centered_z);
  
  // Initialize an EnzoArray representing a component of the interface B-fields.
  // This function accepts both the grouping of permanent and temporary fields,
  // and in each case yields an array that does not include values for the
  // exterior face of the grid. For example, dim=0, the dimension of the
  // resulting EnzoArray is (mz, my-1, mx)
  void load_interior_bfieldi_field(Block *block, Grouping &grouping,
				   int dim, EnzoArray<enzo_float> &array);
  
};
#endif /* ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP */
