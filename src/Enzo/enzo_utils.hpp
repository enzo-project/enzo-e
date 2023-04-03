// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_utils.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-02-11
/// @brief    Generic utility functions useful throughout the Enzo layer

#ifndef ENZO_UTILS_HPP
#define ENZO_UTILS_HPP

/// Namespace for global constants and functions
namespace enzo_utils {

  /// computes the squared magnitude of a 3D vector
  ///
  /// @note
  /// This function uses parentheses to dictate the order of operations.
  /// This is primarily done to take advantage of the `-fprotect-parens` flag
  /// provided by several mainstream c++ compilers (e.g. gcc, clang, icpc).
  /// This flag will honor the order of operations specified by parentheses
  /// even when value-unsafe optimizations for floating-point operations are
  /// enabled (e.g. -ffast-math).
  ///
  /// @note
  /// The choice to annotate this function with FORCE_INLINE was made while
  /// optimizing the Riemann Solvers
  FORCE_INLINE enzo_float squared_mag_vec3D(enzo_float i,
                                            enzo_float j,
                                            enzo_float k) noexcept
  { return ((i*i) + ((j*j) + (k*k))); }

  //----------------------------------------------------------------------

  /// Applies the floor to a quantity.
  ///
  /// This function has primarily been factored out to allow for identifying
  /// erroneous code that causes a floor to be applied in test cases where it
  /// should not be necessary to apply a floor. (This error usually provides a
  /// good starting point for subsequent investigation)
  ///
  /// @note
  /// Maybe this should be a static method of EnzoFluidFloorConfig?
  inline enzo_float apply_floor(const enzo_float value,
                                const enzo_float floor){
    enzo_float out;
    #ifdef RAISE_FLOOR_ERROR
    ASSERT("enzo_utils::apply_floor",
           "DEBUG-MODE ERROR: applying a floor.",
	   value >= floor);
    out = value;
    #else
    out = std::max(value,floor);
    #endif
    return out;
  }

}

#endif /* ENZO_UTILS_HPP */
