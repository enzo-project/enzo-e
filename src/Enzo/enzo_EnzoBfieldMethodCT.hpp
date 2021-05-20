// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBfieldMethodCT.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon June 24 2019
/// @brief    [\ref Enzo] Declaration of BfieldMethodCT

#ifndef ENZO_ENZO_BFIELDMETHODCT_HPP
#define ENZO_ENZO_BFIELDMETHODCT_HPP
class EnzoBfieldMethodCT : public EnzoBfieldMethod
{
  /// @class    EnzoBfieldMethodCT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of constrained transport

public: // interface

  /// Create a new EnzoBfieldMethodCT
  ///
  /// @param[in] num_partial_timesteps The number of partial timesteps over
  ///     which the constructed instance will update the magnetic fields.
  EnzoBfieldMethodCT(int num_partial_timesteps)
    : EnzoBfieldMethod(num_partial_timesteps),
      cell_widths_(nullptr)
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoBfieldMethodCT);

  /// Charm++ PUP::able migration constructor
  EnzoBfieldMethodCT (CkMigrateMessage *m)
    : EnzoBfieldMethod (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p){
    EnzoBfieldMethodCT::pup(p);
  }

  /// checks that the interface bfields exist and have the required shapes
  void check_required_fields() const noexcept;

  /// Overwrites the component of the reconstructed bfields along the axis of
  /// reconstruction with given dimension with corresponding face-centerd
  /// bfield values (which is tracked internally).
  ///
  /// @param[in,out] l_map, r_map These maps should hold the left/right
  ///     reconstructed face-centered magnetic fields. These arrays hold
  ///     face-centered values excluding the exterior faces. The components
  ///     should be associated with the "bfield_x", "bfield_y", and "bfield_z"
  ///     keys.
  /// @param[in]     dim The dimension along which the values were
  ///     reconstructed.
  /// @param[in]     stale_depth The current staling depth. This is the stale
  ///     depth from just before reconstruction plus the reconstructor's
  ///     immediate staling rate.
  void correct_reconstructed_bfield(EnzoEFltArrayMap &l_map,
                                    EnzoEFltArrayMap &r_map, int dim,
                                    int stale_depth) noexcept;

  /// identifies and stores the upwind direction
  ///
  /// @param[in] flux_map Map holding the that holds the density flux along the
  ///     specified dimension in at the "density" key.
  /// @param[in] dim The dimension to identify the upwind direction along.
  /// @param[in] stale_depth The current staling depth. This should match the
  ///     staling depth used to compute the flux_group.
  void identify_upwind(const EnzoEFltArrayMap &flux_map, int dim,
                       int stale_depth) noexcept;

  /// Updates all components of the face-centered and the cell-centered bfields
  ///
  /// @param[in]  cur_prim_map Map containing the current values of the
  ///     cell-centered integrable quantities (before they have been updated
  ///     over the current timestep). Specifically, the velocity and bfield
  ///     components stored in this mapping are used to compute the
  ///     cell-centered E-field.
  /// @param[in]  xflux_map,yflux_map,zflux_map Maps containing the values of
  ///     the fluxes computed along the x, y, and z directions. The function
  ///     namely makes use of the various magnetic field fluxes
  /// @param[out] out_centered_bfield_map Map holding the arrays where the
  ///     updated values for each cell-centered magnetic field component should
  ///     be stored. These arrays can be aliases of `cur_prim_map` (The
  ///     arguments can even reference the same object).
  /// @param[in]  dt The (partial) time-step over which to update the magnetic
  ///     fields.
  /// @param[in]  stale_depth indicates the current stale_depth for the
  ///     supplied quantities. This should nominally be the same as the stale
  ///     depth used to compute the fluxes and that is passed to
  ///     EnzoIntegrableUpdate::update_quantities.
  void update_all_bfield_components(EnzoEFltArrayMap &cur_prim_map,
                                    EnzoEFltArrayMap &xflux_map,
                                    EnzoEFltArrayMap &yflux_map,
                                    EnzoEFltArrayMap &zflux_map,
                                    EnzoEFltArrayMap &out_centered_bfield_map,
                                    enzo_float dt, int stale_depth) noexcept;

  /// Computes a component of the cell-centered magnetic by averaging the
  /// face-centered values of that component of the magnetic field.
  ///
  /// @param[in]  dim The component of the cell-centered B-field to compute.
  ///     Values of 0, 1 and 2 correspond to the x, y and z directions.
  /// @param[out] bfieldc_comp Array where the calculated cell centered values
  ///     of the magnetic field are written.
  /// @param[in]  bfieldi_comp Array containing the interface B-field for the
  ///     `dim` component that are used to compute the cell-centered values.
  ///     This must have the same shape as `bfieldc_comp`, except along
  ///     dimension `dim`. Along that dimension, this array must include one
  ///     more value.
  /// @param[in]  stale_depth indicates the current stale_depth for the
  ///     supplied quantities.
  ///
  /// @note this function is called in `update_all_bfield_components`
  static void compute_center_bfield(int dim, EFlt3DArray &bfieldc_comp,
				    EFlt3DArray &bfieldi_comp,
                                    int stale_depth = 0);

protected: // methods

  /// Virtual method that subclasses overide to preload any method specific
  /// data from the new target block and optionally initialize scratch space
  ///
  /// @param[in] target_block holds data to be processed
  /// @param[in] first_initialization Indicates if this is the first call since
  ///    construction (including deserialization). When true, scratch space
  ///    data should be allocated.
  void register_target_block_(Block *target_block,
                              bool first_initialization) noexcept;

  /// Computes component i of the cell-centered E-field.
  ///
  /// @param[in]  dim The component of the cell-centered E-field to compute.
  ///     Values of 0, 1 and 2 correspond to the x, y and z directions.
  /// @param[out] center_efield The array where the computed cell-centered
  ///     values of the E-field are written.
  /// @param[in]  prim_map Map containing the current values of the
  ///     cell-centered integrable quantities. Specifically, the velocity and
  ///     bfield entries are used to compute the cell-centered E-field.
  /// @param[in]  stale_depth the stale depth at the time of this function call
  ///
  /// @note this function is called in `compute_all_edge_efields`
  static void compute_center_efield(int dim, EFlt3DArray &center_efield,
                                    const EnzoEFltArrayMap &prim_map,
                                    int stale_depth = 0);

  /// Computes component i of the edge-centered E-field that sits on the faces
  /// of dimensions j and k (i, j, and k are any cyclic permutation of x, y, z).
  /// This uses the component i of the cell-centered E-field and component i
  /// of the E-field on the j and k faces (this comes from fluxes computed for
  /// B-field stored in jflux_group and kflux_group). Additionally, it requires
  /// knowledge of the upwind direction on the j and k faces.
  ///
  /// @param[in]  dim The component of the edge-centered E-field to compute.
  ///     Values of 0, 1 and 2 correspond to the x, y and z directions.
  /// @param[in]  center_efield The array holding the cell-centered values for
  ///     component `dim` of the E-field.
  /// @param[out] edge_efield The array where the edge centered values for
  ///     component `dim` of the E-field are to be written. The array should
  ///     have the same number of elements as `center_efield` along dimensions
  ///     `dim` and one more element along the other dimensions (e.g. The
  ///     x-component is cell-centered along the x-direction an face-centered
  ///     along the y- and z-axes).
  /// @param[in]  jflux_map,kflux_map Maps containing the values of the fluxes
  ///     along the j- and k- dimensions (where the dimension i is aligned with
  ///     `dim`). The function namely makes use of the magnetic field fluxes.
  /// @param[in]  weight_l Set of arrays that hold "weight" values. Entries 0,
  ///     1, and 2 are used to hold the x, y, and z components. For a given
  ///     dimension, the array tracks the upwind/downwind direction on the cell
  ///     interfaces for that dimension, without including values on the
  ///     exterior faces of the block (This is included to optionally implement
  ///     the weighting scheme used by Athena++ at a later date).
  /// @param[in]  stale_depth the stale depth at the time of this function call
  ///
  /// @note this function is called in compute_all_edge_efields
  void static compute_edge_efield (int dim, EFlt3DArray &center_efield,
				   EFlt3DArray &edge_efield,
                                   EnzoEFltArrayMap &jflux_map,
                                   EnzoEFltArrayMap &kflux_map,
				   std::array<EFlt3DArray,3> &weight_l,
				   int stale_depth);

  /// Compute the all of the edge-centered electric fields using the current
  /// fluxes and current cell-centered integrable quantities .
  ///
  /// @param[in]  cur_prim_map Map containing the current values of the
  ///     cell-centered integrable quantities (before they have been updated
  ///     over the current timestep). Specifically, the velocity and bfield
  ///     components stored in this mapping are used to compute the
  ///     cell-centered E-field.
  /// @param[in]  xflux_map,yflux_map,zflux_map Maps containing the values of
  ///     the fluxes computed along the x, y, and z directions. The function
  ///     namely makes use of the various magnetic field fluxes
  /// @param[in] center_efield The array where the computed cell-centered
  ///     values of different components of the E-field are temporarily written
  ///     (and then possibly overwritten). This is effectively a preallocated
  ///     scratch buffer whose input values are insignificant.
  /// @param[out] edge_efield_l A set of arrays where the values for each
  ///     component of the edge-centered E-field are to be written. Each entry
  ///     corresponds to a different component. Entries 0, 1, and 2 are used to
  ///     hold the x, y, and z components. If a cell-centered array has shape
  ///     (mz,my,mx), then the arrays should have shapes (mz-1,my-1,mx),
  ///     (mz-1,my,mx-1) and (mz,my-1,mx-1), respectively.
  /// @param[in]  weight_l Set of arrays that hold "weight" values. Entries 0,
  ///     1, and 2 are used to hold the x, y, and z components. For a given
  ///     dimension, the array tracks the upwind/downwind direction on the cell
  ///     interfaces for that dimension, without including values on the
  ///     exterior faces of the block (This is included to optionally implement
  ///     the weighting scheme used by Athena++ at a later date).
  /// @param[in] stale_depth the stale depth at the time of this function call
  static void compute_all_edge_efields
  (EnzoEFltArrayMap &prim_map, EnzoEFltArrayMap &xflux_map,
   EnzoEFltArrayMap &yflux_map, EnzoEFltArrayMap &zflux_map,
   EFlt3DArray &center_efield, std::array<EFlt3DArray,3> &edge_efield_l,
   std::array<EFlt3DArray,3> &weight_l, int stale_depth);

  /// Updates the face-centered B-field component along the ith dimension using
  /// the jth and kth components of the edge-centered E-field
  ///
  /// @param[in]  dim The component of the interface B-field to update.
  ///     Values of 0, 1 and 2 correspond to the x, y and z directions.
  /// @param[in]  edge_efield_l A set of arrays holding the values for each
  ///     component of the edge-centered E-field. Each entry corresponds to a
  ///     different component. Entries 0, 1, and 2 are used to hold the x, y,
  ///     and z components. If a cell-centered array has shape (mz,my,mx), then
  ///     the arrays should have shapes (mz-1,my-1,mx), (mz-1,my,mx-1) and
  ///     (mz,my-1,mx-1), respectively.
  /// @param[in] cur_interface_bfield Array containing the current value of the
  ///     `dim` component of the interface B-field. These values are the values
  ///     from the start of the timestep that are used to compute the updated
  ///     values. Along dimension `dim`, this array should hold one more entry
  ///     than a cell-centered array (it contains values on all faces,
  ///     including the external cell faces). Along the other dimensions, this
  ///     should have the same number of entries as a cell-centered array (e.g.
  ///     the y-component should be cell-centered along x and z, but
  ///     face-centered along y).
  /// @param[out] out_interface_bfield Array where the updated interface value
  ///     of the `dim` component of the interface B-field is written. This
  ///     should have the same shape as `cur_interface_bfield`. This can either
  ///     be the same array passed to `cur_interface_bfield`, or a completely
  ///     different array.
  /// @param[in] dt The time time-step over which to apply the fluxes
  /// @param[in] stale_depth indicates the current stale_depth for the supplied
  ///     quantities
  static void update_bfield(const enzo_float* &cell_widths, int dim,
                            const std::array<EFlt3DArray,3> &efield_l,
                            EFlt3DArray &cur_interface_bfield,
                            EFlt3DArray &out_interface_bfield,
			    enzo_float dt, int stale_depth);

protected: // attributes

  // Block Specific data: (its updated everytime a new block is registered.)

  /// Set of arrays wrapping the Cello fields that contain the interface
  /// magnetic field components. A given component is face-centered along the
  /// face-centered along the dimension of the component. These arrays include
  /// values on the exterior faces of the block.
  ///
  /// If a cell-centered array has shape (mz,my,mx), the entries in this list
  /// have shapes (mz,my,mx+1), (mz,my+1,mx), and (mz+1,my,mx), respectively.
  ///
  /// The entries in this array will be 
  std::array<EFlt3DArray,3> bfieldi_l_;

  /// Array holding cell widths (from the CellWidth attribute of EnzoBlock)
  const enzo_float *cell_widths_;


  // Scratch arrays: The following are reused for different blocks

  /// Set of arrays to temporarily hold the values of the face-centered fields
  /// computed at the half time-step. Each entry should have the same shape as
  /// the corresponding entry in bfieldi_l_
  std::array<EFlt3DArray,3> temp_bfieldi_l_;

  /// Set of arrays that hold "weight" values. For a given dimension, the array
  /// tracks the upwind/downwind direction on the cell interfaces for that
  /// dimension. The arrays exclude values on the exterior faces of the block.
  ///
  /// If a cell-centered array has shape (mz,my,mx), the entries in this list
  /// have shapes (mz,my,mx-1), (mz,my-1,mx), and (mz-1,my,mx), respectively.
  std::array<EFlt3DArray,3> weight_l_;

  /// Set of arrays used to temporarily store the calculated edge-centered
  /// components of the electric field. For a given electric field component,
  /// the corresponding array is cell-centered along that component's dimension
  /// and face-centered along the other dimensions. (e.g. The x-component is
  /// cell-centered along the x-direction an face-centered along the y- and z-
  /// axes)
  ///
  /// If a cell-centered array has shape (mz,my,mx), the entries in this list
  /// have shapes (mz-1,my-1,mx), (mz-1,my,mx-1), and (mz,my-1,mx-1),
  /// respectively.
  std::array<EFlt3DArray,3> edge_efield_l_;

  /// Array used to temporarily store cell-centered values of different
  /// electric field components. This array is reused for each component.
  EFlt3DArray center_efield_;

};
#endif /* ENZO_ENZO_BFIELDMETHODCT_HPP */
