// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFieldArrayFactory.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implementation of Enzo's FieldArrayFactory. It
///           streamlines the loading of fields as arrays

#ifndef ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP
#define ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP

class EnzoFieldArrayFactory
{
  /// @class    EnzoFieldArrayFactory
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates construction of CelloArrays from fields

public:

  /// Create a new EnzoFieldArrayFactory
  EnzoFieldArrayFactory(Block *block, int stale_depth = 0)
  {
    block_ = block;
    ASSERT("EnzoFieldArrayFactory", "stale_depth must be >= 0",
	   stale_depth >= 0);
    stale_depth_ = stale_depth;
  }

  ~EnzoFieldArrayFactory()
  { block_ = NULL;}

  /// returns a field as an array
  EFlt3DArray from_name(std::string field_name);

  /// returns a field from a Grouping
  EFlt3DArray from_grouping(Grouping &grouping, std::string group_name,
			    int index);

  /// Read in fields from Grouping that represented reconstructed quantities.
  /// The reconstructed fields are formally cell-centered fields. However, in
  /// reality they act as face-centered fields (that exclude faces on the
  /// exterior of the grid) that get reused for different dimensions. This
  /// function constructs the array such that it is face-centered along
  /// dimension dim
  EFlt3DArray reconstructed_field(Grouping &grouping, std::string group_name,
				  int index, int dim);

  /// Read in field from Grouping of face-centered interface B-fields. The
  /// returned view doesn't include face-centered values on the exterior of
  /// the grid.
  EFlt3DArray interior_bfieldi(Grouping &grouping, int dim);
protected: // methods

  /// Helper function that checks the validity of a group_name and index
  void check_grouping_details_(Grouping &grouping, std::string group_name,
			       int index);
  
  /// Helper function that returns an instance of EFlt3DArray but without the
  /// stale cells (stale_depth must be a positive integer)
  EFlt3DArray exclude_stale_cells_(EFlt3DArray &arr);

protected: // attributes

  /// Contains the relevant current data
  Block* block_;

  /// indicates the number of field entries from the outermost field value that
  /// the region including "stale" values extends over.
  int stale_depth_;
};
#endif /* ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP */
