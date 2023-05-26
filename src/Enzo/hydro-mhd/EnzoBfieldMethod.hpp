// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBfieldMethod.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Tues April 20 2021
/// @brief    [\ref Enzo] Declaration of the BfieldMethod abstract base class.
///           This class should be subclassed to implement different magnetic
///           field integration methods.

#ifndef ENZO_ENZO_BFIELDMETHOD_HPP
#define ENZO_ENZO_BFIELDMETHOD_HPP
class EnzoBfieldMethod : public PUP::able
{
  /// @class    EnzoBfieldMethod
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates methods for integrating magnetic
  ///           fields.
  ///
  /// Subclass instances are meant to be used as components of an MHD or
  /// hydrodynamical integrators that are entirely responsible for
  /// bfield integration-related operations that would not otherwise be
  /// performed by the solver (e.g. constrained transport, divergence cleaning)
  /// Note that integration-related operations that are similar to ordinary
  /// hydrodynamic integration (e.g. reconstruction and flux calculation) will
  /// still be performed by the integrator. The distinction between whether a
  /// pointer is a `nullptr` or points to an instance of this class should
  /// entirely dictate whether these extra integration operations occur.
  ///
  /// To accomplish these goals, subclass instances are responsible for the
  /// allocation/deallocation of the required temporary scratch-space
  /// and providing the actual integration-related operations
  ///
  /// This class effectively implements a state machine.
  /// - To update the bfield for a given block, the `register_target_block`
  ///   method must first be used to register that block. An error will be
  ///   raised if this method is called and a block has already been
  ///   initialized. Because the shapes of the blocks are not necessarily known
  ///   at construction, the scratch arrays are lazily initialized the first
  ///   time this is called.
  /// - Errors will be raised if the `correct_reconstructed_bfield`,
  ///   `identify_upwind`, `update_all_bfield_components`, or
  ///   `increment_partial_timestep` are called without having a registered
  ///    block.
  /// - When the `increment_partial_timestep` method has been called
  ///   `num_partial_timesteps` times after a block has been registered, the
  ///   block will be unregistered.
  ///
  /// To the hydro solver the main effect observable effect of calling this
  /// class is updating the cell-centered Bfield components. However, the side
  /// effects of the operations implicity performed by this class may be
  /// equally as important. For example, the subclass for Constrained Transport 
  /// needs to updates the face-centered Bfield components (given that the
  /// face-centered fields are the primary represntation of the Bfield and the
  /// cell-centered values are derived directly from the face-centered values).
  

public: // interface

  /// Create a new EnzoBfieldMethod object.
  ///
  /// @param[in] num_partial_timesteps The number of partial timesteps over
  ///     which the constructed instance will update the magnetic fields.
  EnzoBfieldMethod(int num_partial_timesteps);

  /// Virtual destructor
  virtual ~EnzoBfieldMethod()
  {}

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoBfieldMethod);

  /// CHARM++ migration constructor for PUP::able
  EnzoBfieldMethod (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  ///
  /// Subclasses shouldn't PUP scratch arrays.
  void pup (PUP::er &p);

  // State machine Methods:

  /// Register the pointer to the Block for which the Bfields will be updated.
  ///
  /// Only one block can be registered at a time (an error will be raised if
  /// another block is already registered). The very first time this is called,
  /// it will initialize the scratch space.
  ///
  /// @param[in] block holds data to be processed
  void register_target_block(Block *block) noexcept;

  /// Returns the partial timestep index.
  ///
  /// 0 means that the current interface bfields are stored in the permanent
  /// fields. 1 means that they are stored in the temporary arrays. An error
  /// will occur if no block is registered.
  int partial_timestep_index() const noexcept
  {
    require_registered_block_();
    return partial_timestep_index_;
  }

  /// increments the partial timestep index
  void increment_partial_timestep() noexcept;


  // Descriptor methods

  /// checks that all required fields exist and have the required shapes
  virtual void check_required_fields() const noexcept = 0;

  // Physics methods
  
  /// In the case of Constrained Transport, overwrites the component of the
  /// reconstructed bfield corresponding to `dim` with the corresponding
  /// face-centered bfield values (which is tracked internally).
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
  virtual void correct_reconstructed_bfield(EnzoEFltArrayMap &l_map,
                                            EnzoEFltArrayMap &r_map, int dim,
                                            int stale_depth) noexcept = 0;

  /// In the case of Constrained Transport, identifies and stores the upwind
  /// direction.
  ///
  /// @param[in] flux_map Map holding the that holds the density flux along the
  ///     specified dimension in at the "density" key.
  /// @param[in] dim The dimension to identify the upwind direction along.
  /// @param[in] stale_depth The current staling depth. This should match the
  ///     staling depth used to compute the flux_group.
  virtual void identify_upwind(const EnzoEFltArrayMap &flux_map, int dim,
                               int stale_depth) noexcept = 0;

  /// Updates all components of the bfields (this is to be called before the
  /// hydro quantities are updated)
  ///
  /// For constrained transport, this includes both the face-centered and the
  /// cell-centered bfields (the face-centered bfield updates are handled
  /// completely internally). The velocity and bfield components stored in this
  /// `cur_integration_map` map are used to compute the cell-centered E-field.
  ///
  /// @param[in]  cur_integration_map Map containing the current values of the
  ///     cell-centered integration quantities (before they have been updated
  ///     over the current timestep). 
  /// @param[in]  xflux_map,yflux_map,zflux_map Maps containing the values of
  ///     the fluxes computed along the x, y, and z directions. The function
  ///     namely makes use of the various magnetic field fluxes
  /// @param[out] out_centered_bfield_map Map holding the arrays where the
  ///     updated values for each cell-centered magnetic field component should
  ///     be stored. This can be an alias of `cur_integration_map`.
  /// @param[in]  dt The (partial) time-step over which to update the magnetic
  ///     fields.
  /// @param[in]  stale_depth indicates the current stale_depth for the
  ///     supplied quantities. This should nominally be the same as the stale
  ///     depth used to compute the fluxes and that is passed to
  ///     EnzoIntegrationQuanUpdate::update_quantities.
  virtual void update_all_bfield_components
  (EnzoEFltArrayMap &cur_integration_map, const EnzoEFltArrayMap &xflux_map,
   const EnzoEFltArrayMap &yflux_map, const EnzoEFltArrayMap &zflux_map,
   EnzoEFltArrayMap &out_centered_bfield_map, enzo_float dt, int stale_depth)
    noexcept=0;

protected:

  /// Virtual method that subclasses overide to preload any method specific
  /// data from the new target block and optionally initialize scratch space
  ///
  /// @param[in] target_block holds data to be processed
  /// @param[in] first_initialization Indicates if this is the first call since
  ///    construction (including deserialization). When true, scratch space
  ///    data should be allocated.
  virtual void register_target_block_(Block *target_block,
                                      bool first_initialization) noexcept = 0;

  /// Helper method that raises an error if there isn't a registered block.
  void require_registered_block_() const noexcept
  {
    if (target_block_ == nullptr){
      ERROR("EnzoConstrainedTransport::require_registered_block_",
            "An invalid method has been executed that requires that a target "
            "block is registered.");
    }
  }

  int num_partial_timesteps() const noexcept {return num_partial_timesteps_;}

private: // keep subclasses from touching these attributes

  /// Number of partial timesteps in a cycle. A value of 1 would means that
  /// there is only a single update over the full timestep. A value of 2 means
  /// that there is a half timestep and a full timestep.
  int num_partial_timesteps_;

  /// Holds the index of the current partial timestep. This will be reset
  /// everytime a new target_block_ is registered. A negative value is used
  /// to indicate that the scratch arrays are not allocated.
  int partial_timestep_index_;

  /// Pointer to the Block that this class is currently operating upon. To try
  /// to keep the physics as decoupled as possible, we avoid letting the
  /// subclasses from directly touching the target_block_
  Block *target_block_;

};

#endif /* ENZO_ENZO_BFIELDMETHOD_HPP */

