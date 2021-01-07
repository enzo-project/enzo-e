// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoTempFieldUtils.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Tues April 28 2020
/// @brief    [\ref Enzo] Declaration of the EnzoTempFieldUtils namespace

#ifndef ENZO_TEMP_FIELDS_UTILS_HPP
#define ENZO_TEMP_FIELDS_UTILS_HPP

namespace EnzoTempFieldUtils {
  /// @namespace EnzoTempFieldUtils
  /// @ingroup   Enzo
  /// @brief     [\ref Enzo] Collection of functions that are useful for
  ///            repeated creation and deletion of temporary fields.
  ///
  /// This solution is somewhat inelegant. Ideally in the future, public use of
  /// these functions will be replaced with a Context Manager Class whose
  /// lifetime is coupled to the allocation and deallocation of reused
  /// temporary fields (OR local temporary fields are not used).

  /// Allocates a frequently reused temporary field
  ///
  /// If the field was not previously created it is inserted into
  /// field.field_descr. If the field_name matches the name of the permanent
  /// field, then nothing happens.
  ///
  /// @param field Field object that will hold the allocated temporary field
  /// @param field_name The name of the temporary field to allocate
  /// @param cx,cy,cz The centering of the temporary field to allocate
  void prep_reused_temp_field(Field &field, std::string field_name,
			      int cx, int cy, int cz);

  /// deallocates the temporary fields listed in a grouping
  ///
  /// @param field Field object that holds the temporary fields to be
  ///     deallocated
  /// @param group_names Vector of group names that indicates which groups the
  ///     fields that should be deallocated are a part of.
  /// @param grouping Together with group_names, this specified which temporary
  ///     fields should be deallocated.
  void deallocate_grouping_fields(Field &field,
				  const std::vector<std::string> &group_names,
				  Grouping &grouping);

}
#endif /* ENZO_TEMP_FIELDS_UTILS_HPP */
