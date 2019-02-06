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

  // read in Grouping fields
  EFlt3DArray from_grouping(Grouping &grouping, std::string group_name,
			    int index);

  // read in fields from Grouping represented reconstructed quantities
  //
  // Reconstructed fields are cell-centered fields that do not include faces on
  // the exterior of the grid. Because the reconstructed fields are reused for
  // different dimensions, they are not registerred internally as being
  // face-centered. The dim argument indicates the direction along which the
  // field should be face-centered.
  EFlt3DArray reconstructed_field(Grouping &grouping, std::string group_name,
				  int index, int dim);

  // read in field from Grouping of face-centered interface B-fields. The
  // resulting view does not include face-centered values on the exterior of
  // the grid.
  EFlt3DArray interior_bfieldi(Grouping &grouping, int dim);

  // Depreciated.
  //
  // The interface was written this way due to a misconception.
  // Most recently this was used for debugging and for reconstructed fields.
  //
  // Setting cell_centered_{dim} to true indicates that values
  // are cell-centered along {dim}
  EFlt3DArray load_temp_interface_grouping_field(Grouping &grouping,
						 std::string group_name,
						 int index,
						 bool cell_centered_x,
						 bool cell_centered_y,
						 bool cell_centered_z);

private:
  Block* block_;
};
#endif /* ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP */
