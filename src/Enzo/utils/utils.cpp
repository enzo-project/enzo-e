// See LICENSE_CELLO file for license and copyright information

/// @file     utils.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2024-02-20
/// @brief    Definitions of generic utility functions useful throughout the
///           Enzo layer

#include "Enzo/enzo.hpp"
#include "Enzo/utils/utils.hpp"

#include <cmath>

//----------------------------------------------------------------------

static enzo_float median(enzo_float x, enzo_float y, enzo_float z) noexcept
{
  // we intentionally use bitwise-and rather than logical-and for speed
  if ((x >= y) & (x <= z)) {
    return x;
  } else if ((x >= z) & (x <= y)) {
    return x;
  } else if ((y >= x) & (y <= z)) {
    return y;
  } else if ((y >= z) & (y <= x)) {
    return y;
  }
  return z;
}

bool enzo_utils::consistent_cube_cellwidths(enzo_float dx, enzo_float dy,
                                            enzo_float dz) noexcept
{
  // we want to ensure the median value is within N ULPs of the other 2 values
  // - a ulp is the difference between 2 adjacent floating point values
  //
  // We'll start with N = 4. But we may want need to change that

  enzo_float width = median(dx,dy,dz);

  auto helper = [width](enzo_float toward) -> enzo_float
  {
    enzo_float cur = width;
    for (int i = 0; i < 4; i++){ cur = std::nextafter(cur, toward); }
  };

  enzo_float min_val = helper(-INFINITY);
  enzo_float max_val = helper(+INFINITY);

  // we intentionally use bitwise-and rather than logical-and for speed
  return ((min_val <= dx) & (min_val <= dy) & (min_val <= dz) &
          (dx <= max_val) & (dy <= max_val) & (dz <= max_val));
}
