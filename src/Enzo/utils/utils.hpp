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

  /// compute the minimum of 3 values
  ///
  /// @note
  /// adapted from Enzo's ReconstructionRoutines.h
  template <typename T>
  T min(T a, T b, T c){
    if (a<b) {
      return (c<a) ? c : a;
    } else {
      return (c<b) ? c : b;
    }
  }

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

  //----------------------------------------------------------------------

  /// Utiltity template function used for executing a Compute-Kernel (currently
  /// only used for CPU bound tasks)
  ///
  /// This is most useful when K is a lambda
  template<class K>
  void exec_loop(int mz, int my, int mx, int stale_depth, K& kernel) {
    const int rank = cello::rank();

    const int ix_start = stale_depth;
    const int ix_stop = mx - stale_depth;

    const int iy_start = (rank > 1) ? stale_depth : 0;
    const int iy_stop  = (rank > 1) ? my - stale_depth : my;

    const int iz_start = (rank > 2) ? stale_depth : 0;
    const int iz_stop  = (rank > 2) ? mz - stale_depth : mz;

    if ((ix_start >= ix_stop) | (iy_start >= iy_stop) | (iz_start >= iz_stop)){
      ERROR("enzo_utils::exec_loop", "stale_depth is too large");
    } else if (stale_depth < 0){
      ERROR("enzo_utils::exec_loop", "stale_depth is negative");
    }

    for (int iz = iz_start; iz < iz_stop; iz++){
      for (int iy = iy_start; iy < iy_stop; iy++){
        for (int ix = ix_start; ix < ix_stop; ix++){
          kernel(iz, iy, ix);
        }
      }
    }
  }

  //----------------------------------------------------------------------

  /// checks whether the existence of certain velocity fields are consistent
  /// with the rank. If it isn't consistent, the program aborts with an
  /// informative error
  ///
  /// This function allows the ``"velocity_x"``, ``"velocity_y"``, and
  /// ``"velocity_z"`` to exist when ``cello::rank()`` is 1, 2, or 3. It also
  /// permits the following 2 cases:
  ///  * when ``cello::rank() == 1``, just the ``"velocity_x"`` field exists 
  ///  * when ``cello::rank() == 2``, just the ``"velocity_x"`` &
  ///   ``"velocity_y"`` fields exist.
  void assert_rank_velocity_field_consistency(FieldDescr&) noexcept;

}

#endif /* ENZO_UTILS_HPP */
