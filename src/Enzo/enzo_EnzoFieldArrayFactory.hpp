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
  /// @brief    [\ref Enzo] Encapsulates construction of CelloArrays that wrap
  ///           Cello fields
  ///
  /// This functionallity should probably be incorporated directly into Cello's
  /// Field class

public:

  /// Create a new EnzoFieldArrayFactory
  EnzoFieldArrayFactory(Block *block, int stale_depth = 0)
  {
    ASSERT("EnzoFieldArrayFactory", "block must not be a nullptr",
	   block != nullptr);
    block_ = block;
    ASSERT("EnzoFieldArrayFactory", "stale_depth must be >= 0",
	   stale_depth >= 0);
    stale_depth_ = stale_depth;
  }

  /// returns a field as an array
  EFlt3DArray from_name(const std::string &field_name){
    Field field = block_->data()->field();
    const int id = field.field_id(field_name);

    ASSERT1("EnzoFieldArrayFactory",
            "The \"%s\" field has the incorrect precision.",
            field_name.c_str(),
            cello::sizeof_precision(field.precision(id)) == sizeof(enzo_float));

    int mx, my, mz;
    field.dimensions (id,&mx,&my,&mz);
    EFlt3DArray arr((enzo_float *) field.values(field_name), mz, my, mx);

    if (stale_depth_ != 0){
      ASSERT("EnzoFieldArrayFactory",
             "each dim of arr must exceed 2*stale_depth_.",
             2*stale_depth_ < arr.shape(0) &&
             2*stale_depth_ < arr.shape(1) &&
             2*stale_depth_ < arr.shape(2));
    return arr.subarray(CSlice(stale_depth_, -stale_depth_),
                        CSlice(stale_depth_, -stale_depth_),
                        CSlice(stale_depth_, -stale_depth_));
    }
    return arr;
  }

protected: // attributes

  /// Contains the relevant current data
  Block* block_;

  /// indicates the number of field entries from the outermost field value that
  /// the region including "stale" values extends over.
  int stale_depth_;
};
#endif /* ENZO_ENZO_FIELD_ARRAY_FACTORY_HPP */
