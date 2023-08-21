// See LICENSE_CELLO file for license and copyright information

/// @file     utils.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-08-21
/// @brief    Implements generic utility functions useful throughout Enzo layer

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void enzo_utils::assert_rank_velocity_field_consistency(FieldDescr& f_descr)
  noexcept
{
  const int rank = cello::rank();

  if ((rank < 1) || (rank > 3)) {
    ERROR("enzo_utils::assert_rank_velocity_field_consistency",
          "something is horribly wrong: cello::rank() must return 1, 2, or 3");
  }

  const bool has_vx = f_descr.is_field("velocity_x");
  const bool has_vy = f_descr.is_field("velocity_y");
  const bool has_vz = f_descr.is_field("velocity_z");

  if (!has_vx) {
    ERROR("enzo_utils::assert_rank_velocity_field_consistency",
          "\"velocity_x\" field must be defined for any rank");
  } else if ((rank == 1) && (has_vy != has_vz)) {
    ERROR("enzo_utils::assert_rank_velocity_field_consistency",
          "when rank == 1, the \"velocity_y\" and \"velocity_z\" fields can "
          "both be defined. Alternatively, neither should be defined.");
  } else if ((!has_vy) && (rank > 1)) {
    ERROR("enzo_utils::assert_rank_velocity_field_consistency",
          "\"velocity_y\" field must be defined when rank > 1");
  } else if ((!has_vy) && has_vz) {
    ERROR("enzo_utils::assert_rank_velocity_field_consistency",
          "velocity_y field must be defined velocity_z field is defined");
  } else if ((!has_vz) && rank > 2) {
    ERROR("enzo_utils::assert_rank_velocity_field_consistency",
          "velocity_z field must be defined when rank > 2");
  }
}
