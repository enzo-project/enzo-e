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
  // we want to ensure the median of (dx,dy,dz) is within N ULPs of the other 2
  // values
  // - A ULP, or 'Unit in the Last Place', is the difference between 2 adjacent
  //   floating point values.
  // - Since floating point calculations are not exact and decimal numbers
  //   cannot always be represented exactly, comparisons between two nominally
  //   identical numbers need to take that into account.
  // - To that end, this routine checks to make sure all values are very close
  //   to each other - within N ULPs. Thus it is checking to make sure the
  //   three values are approximately equal instead of exactly equal.
  //
  // We currently define approximate equality to be within 4 ULPs of the median
  // of the three values (as can be seen by the value of N_ULP)
  // - This may need to be changed if this tolerance ends up being too
  //   restrictive.
  const int N_ULP = 4;

  // find the median cell-width
  const enzo_float median_cellwidth = median(dx,dy,dz);

  // in the following chunk of code, we search for 2 values:
  // 1. a value that is `N_ULP` smaller than median_cellwidth
  // 2. a value that is `N_ULP` larger than median_cellwidth
  //
  // To accomplish this:
  // -> we initialize the variables `lower_val` and `upper_val` and initially
  //    assign them the value stored in `median_cellwidth`
  // -> each time we pass through the for-loop, we update `lower_val` so it
  //    stores the next representable floating point that is smaller than its
  //    current value. At the same time, we also update `upper_val` so that it
  //    stores the representable floating point that is smaller than its
  //    current value. (Essentially, we change each variable by 1 ULP each time
  //    we pass through the loop)
  enzo_float lower_val = median_cellwidth;
  enzo_float upper_val = median_cellwidth;
  for (int i = 0; i < N_ULP; i++){
    // we use `-INFINITY` and `+INFINITY` to specify whether std::nextafter
    // should returns a smaller or larger value
    lower_val = std::nextafter(lower_val, -INFINITY);
    upper_val = std::nextafter(upper_val, +INFINITY);
  }

  // For explicitness:
  // - `min_val` is `N_ULP` ULPs smaller than `median_cellwidth`
  // - `max_val` is `N_ULP` ULPs larger than `median_cellwidth`
  const enzo_float min_val = lower_val;
  const enzo_float max_val = upper_val;

  // The function returns true when `dx`, `dy`, and `dz` all have values within
  // the interval of values bounded by `min_val <= val <= upper_val`
  // - in other words, we return `true` if all 3 values are within `N_ULP` ULPs
  //   of each other. Otherwise, we return `false`.
  // - we intentionally use bitwise-and rather than logical-and for speed
  return ((min_val <= dx) & (dx <= max_val) &
          (min_val <= dy) & (dy <= max_val) &
          (min_val <= dz) & (dz <= max_val));
}
