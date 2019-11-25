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

  /// returns a field from a Grouping as an array with a shape appropriate
  /// for being face-centered along the specified axis (excluding faces on the
  /// exterior of the grid).
  ///
  /// If the underlying field is cell-centered, then the returned array wraps
  /// an appropriate subset of the underlying data. This facillitates the reuse
  /// of individual fields for temporarily storing data that is face-centered
  /// along different axes (e.g. like for reconstructed fields)
  ///
  /// If the underlying field is face-centered, then the field is returned as
  /// an array (equivalent to calling `from_grouping`) as long as the field's
  /// properties exactly match the designated centering (i.e. it is only
  /// face-center along the specified axis AND it excludes the faces on the
  /// grid exterior). If the properties don't match an error is raised.  
  EFlt3DArray assigned_center_from_grouping(Grouping &grouping,
					    std::string group_name,
					    int index, int dim);

  /// returns a field as an array with a shape appropriate for being
  /// face-centered along a specified axis.
  EFlt3DArray assigned_center_from_name(std::string field_name, int dim)
  {
    Grouping temp_group;
    temp_group.add(field_name, "group");
    return assigned_center_from_grouping(temp_group, "group", 0, dim);
  }

  /// Read in fields from Grouping that represented reconstructed quantities.
  /// The reconstructed fields should formally be initialized as cell-centered
  /// fields.
  ///
  /// This is deprecated and replaced by `assigned_center_from_grouping`
  EFlt3DArray reconstructed_field(Grouping &grouping, std::string group_name,
				  int index, int dim)
  { return assigned_center_from_grouping(grouping, group_name, index, dim); }

  /// Read in field from Grouping of face-centered interface B-fields. The
  /// returned view doesn't include face-centered values on the exterior of
  /// the grid.
  EFlt3DArray interior_bfieldi(Grouping &grouping, int dim);
protected: // methods

  /// Helper function that reads in the field without applying stale depth
  EFlt3DArray full_field_from_name_(std::string field_name);

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
