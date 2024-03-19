// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructor.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implementation of Enzo's Reconstructor interface

#ifndef ENZO_ENZO_RECONSTRUCTOR_HPP
#define ENZO_ENZO_RECONSTRUCTOR_HPP

class EnzoReconstructor
{
  /// @class    EnzoReconstructor
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates reconstruction of primitives at
  ///           cell interfaces

public: // interface

  /// Factory method for constructing EnzoReconstructor
  /// (The signature of this method may need to be modified)
  ///
  /// @param[in] active_primitive_keys A vector listing the names for each
  ///     primitive quantity that are to be reconstructed. This should omit the
  ///     names of the passively advected scalars.
  /// @param[in] name The name of the Riemann solver to use. Valid names
  ///     include "nn" and "plm"
  /// @param[in] theta_limiter An argument that is optionally used to tune
  ///     certain types of limiters
  static EnzoReconstructor* construct_reconstructor
    (const std::vector<std::string> &active_primitive_keys,
     std::string name, enzo_float theta_limiter);

  /// Create a new EnzoReconstructor
  EnzoReconstructor(std::vector<std::string> active_key_names) throw()
    : active_key_names_(active_key_names)
  { }

  /// Virtual destructor
  virtual ~EnzoReconstructor()
  {  }

  /// Reconstructs the interface values
  ///
  /// @param[in]  prim_map Map holding the data for the cell-centered
  ///     primitives that are to be reconstructed. This is expected to have
  ///     a key-array pair for each entry in the list passed as the
  ///     ``active_primitive_keys`` argument of the factory method
  ///     (plus all of the keys listed in `passive lists`)
  /// @param[out] priml_map,primr_map Holds existing arrays where the
  ///     left/right reconstructed, face-centered primitives are written.
  ///     These must supply the same keys that are expected for prim_map.
  ///     The arrays are expected to have identical shapes to those in
  ///     prim_map, except along dimension `dim`. Along that dimension, these
  ///     arrays should hold one fewer value.
  /// @param[in]  dim Dimension along which to reconstruct interface values.
  ///     Values of 0, 1, and 2 correspond to the x, y, and z directions.
  /// @param[in]  stale_depth indicates the current stale_depth for the
  ///     supplied cell-centered quantities
  /// @param[in]  passive_list A list of keys for passive scalars.
  virtual void reconstruct_interface
  (const EnzoEFltArrayMap &prim_map, EnzoEFltArrayMap &priml_map,
   EnzoEFltArrayMap &primr_map, const int dim, const int stale_depth,
   const str_vec_t& passive_list)=0;

  /// The rate amount by which the stale_depth increases after the current
  /// reconstructor is used to update the fluid over a (partial or full)
  /// time-step. If the fluid is update over a partial timestep before being
  /// updated, increases to the stale_depth are cummulative.
  ///
  /// stale_depth indicates the number of field entries from the outermost
  /// field value that the region including "stale" values (need to be
  /// refreshed) extends over. For cell-centered fields this is the number
  /// of cells away from the edge. Along the dimension(s) of face-centering,
  /// the region of valid (un-staled) entries for a face-centered field,
  /// with/without values on the exterior of the grid, will ALWAYS contain 1
  /// more/less valid field entry than the region of valid entries for a
  /// cell-centered field
  /// 
  /// the staling_rate can be decomposed into 2 parts: immediate and delayed
  /// - immediate staling rate is the amount by which stale_depth increases
  ///   immediately after performing reconstruction.
  /// - delayed staling rate is the amount by which stale_depth increases
  ///   after the fluid is updated over a (partial or full) time-step using the
  ///   fluxes computed from the reconstructed values.
  virtual int total_staling_rate()=0;

  /// immediate staling rate is the amount by which stale_depth increases
  /// immediately after performing reconstruction. This is equal to the number
  /// cells separating the first (last) cell interface that doesn't have both a
  /// left AND right reconstructed values from the edge of the grid. For
  /// interpolation like nearest neighbor, this is 0, while for interpolation
  /// like piecewise linear, this is 1.
  virtual int immediate_staling_rate()=0;

  /// delayed staling rate is the amount by which stale_depth increases after
  /// the fluid is updated over a (partial or full) time-step using the fluxes
  /// computed from the reconstructed values.
  int delayed_staling_rate()
  { return total_staling_rate() - immediate_staling_rate(); }

protected:
  /// list of the key names for all components of (non-passively advected)
  /// primitives that are to be reconstructed.
  std::vector<std::string> active_key_names_;
};

#endif /* ENZO_ENZO_RECONSTRUCTOR_HPP */
