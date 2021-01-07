// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoConstrainedTransport.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon June 24 2019
/// @brief    [\ref Enzo] Implementation of EnzoConstrainedTransport

#ifndef ENZO_ENZO_CONSTRAINEDTRANSPORT_HPP
#define ENZO_ENZO_CONSTRAINEDTRANSPORT_HPP
class EnzoConstrainedTransport
{
  /// @class    EnzoConstrainedTransport
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of constrained transport
  ///
  /// Instances of this class are meant to be used as components of an MHD or
  /// hydrodynamical integrators that are entirely responsible for all
  /// constrained-transport related operations that would not otherwise be
  /// performed by the solver. The distinction between whether a pointer is
  /// NULL or points to an instance of this class should entirely dictate
  /// whether constrained-transport related operations are used in the
  /// integrator.
  ///
  /// To accomplish these goals, this object is responsible for the
  /// allocation/deallocation of the required temporary scratch-space
  /// and providing the actual constrained transport operations. To manage
  /// scratch-space in a thread-safe manner (i.e. allocating/deallocating it
  /// separately for each block) while avoiding two-step constructing and
  /// destruction of instances of this class, the temporary scratch-space
  /// arrays are implemented as temporary fields that are allocated/deallocated
  /// in the class's constructor/destructor. As a consequence, instances of
  /// this class need to be constructed and destroyed separately for each leaf
  /// block during every cycle that the underlying MHD integrator is applied.
  /// The instance also internally tracks the current partial timestep to be
  /// less intrusive for the underlying MHD integrator.
  ///
  /// Although this design is crude, it is simple. The design could be improved
  /// by making separating the responsiblities of the allocation/deallocation
  /// of temporary scratch space and the implementation of the constrained
  /// transport operation between two classes. However this refactoring has
  /// been defered to the future when other classes may be introduced to
  /// encapsulate other approaches to updated the magnetic field (e.g. Dedner
  /// divergence cleaning).

public: // interface

  /// Create a new EnzoConstrainedTransport (allocates scratch space)
  ///
  /// @param block holds data to be processed
  /// @param num_partial_timesteps The number of partial timesteps over which
  ///     the constructed instance will update the magnetic fields.
  EnzoConstrainedTransport(Block *block, int num_partial_timesteps);

  /// Destructor (deallocates the scratch space)
  ~EnzoConstrainedTransport();

  /// adds the interface bfields to the refresh list (and makes sure that they
  /// exist)
  static void update_refresh(Refresh* refresh);

  /// Returns the partial timestep index.
  ///
  /// 0 means that the current interface bfields are stored in the permanent
  /// fields (this is how the. 1 means that they are stored in the 
  int partial_timestep_index() const { return partial_timestep_index_; }

  /// increments the partial timestep index
  void increment_partial_timestep() throw();

  /// Overwrites the component of the reconstructed bfields along the axis of
  /// reconstruction with given dimension with corresponding face-centerd
  /// bfield values (which is tracked internally).
  ///
  /// @param l_group, r_group These should each contain a group called "bfield"
  ///     that holds the names of the fields which store the components of the
  ///     left/right reconstructed face-centered magnetic fields. The relevant
  ///     fields should be formally defined as cell-centered. During the
  ///     calculation, they are treated as face-centered (without having values
  ///     on the exterior faces of the block).
  ///     holding the hold the 
  /// @param dim The dimension along which the values were reconstructed.
  /// @param stale_depth The current staling depth. This is the stale depth
  ///     from just before reconstruction plus the reconstructor's immediate
  ///     staling rate.
  void correct_reconstructed_bfield(Grouping &l_group, Grouping &r_group,
				    int dim, int stale_depth);

  /// identifies and stores the upwind direction
  ///
  /// @param flux_group this must have a "density" group that contains one name
  ///     that refers to the field holding the "density" flux along the
  ///     specified dimemsion.
  /// @param dim The dimension to identify the upwind direction along.
  /// @param stale_depth The current staling depth. This should match the
  ///     staling depth used to compute the flux_group.
  void identify_upwind(Grouping &flux_group, int dim, int stale_depth);

  /// Updates all components of the face-centered and the cell-centered bfields
  ///
  /// @param cur_prim_group Grouping of the fields containing the current
  ///     values of the cell-centered integrable quantities (before dt is
  ///     applied). The cell-centered E-field is computed from the stored
  ///     bfield and velocity fields (which should be stored in the "bfield"
  ///     and "velocity" group).
  /// @param xflux_group,yflux_group,zflux_group holds field names where the
  ///     fluxes along the x, y, and z directions are stored. These should
  ///     include fluxes computed for the magnetic fields
  /// @param bfieldc_group this must have a group called "bfield" holding 3
  ///     cell-centered fields (one for each B-field component). This can be
  ///     the same as cur_prim_group
  /// @param dt The (partial) time-step over which to apply the fluxes
  /// @param stale_depth indicates the current stale_depth for the supplied
  ///     quantities. This should nominally be the same stale depth as the
  ///     value used to compute the fluxes and passed to
  ///     EnzoIntegrableUpdate::update_quantities.
  void update_all_bfield_components(Grouping &cur_prim_group,
				    Grouping &xflux_group,
				    Grouping &yflux_group,
				    Grouping &zflux_group,
				    Grouping &out_centered_bfield_group,
				    enzo_float dt, int stale_depth);

  /// Computes a component of the cell-centered magnetic by averaging the
  /// face-centered values of that component of the magnetic field.
  ///
  /// @param block holds data to be processed
  /// @param dim The component of the cell-centered B-field to compute. Values
  ///     of 0, 1 and 2 correspond to the x, y and z directions, respectively.
  /// @param bfieldc_group this must have a group called "bfield" holding 3
  ///     cell-centered fields (one for each B-field component). The field
  ///     holding the data for the component corresponding to dim will have its
  ///     values updated by this method.
  /// @param bfieldi_group contains the interface B-field used to update the
  ///     cell-centered value. This must have a group called "bfield" holding 3
  ///     face-centered fields (one for each B-field component). A given
  ///     component should be cell-centered for the dimensions that don't match
  ///     the name of the component and face-centered (including space for
  ///     values located on exterior faces of the block) along the dimension
  ///     that matches the name of the component (e.g. the y-component should
  ///     be cell-centered along x and z, but face-centered along y).
  /// @param stale_depth indicates the current stale_depth for the supplied
  ///     quantities
  static void compute_center_bfield(Block *block, int dim,
				    Grouping &bfieldc_group,
				    Grouping &bfieldi_group,
				    int stale_depth = 0);

protected: // methods
  /// Computes component i of the cell-centered E-field.
  ///
  /// @param block holds data to be processed
  /// @param dim The component of the cell-centered E-field to compute. Values
  ///     of 0, 1 and 2 correspond to the x, y and z directions, respectively.
  /// @param center_efield_name The name of the field that will store the
  ///     calculated values of the cell-centered E-field
  /// @param prim_group holds the fields containing the current values of the
  ///     cell-centered integrable quantities (E-field is computed from the
  ///     stored bfield and velocity)
  /// @param stale_depth the stale depth at the time of this function call
  ///
  /// @note this function is called in compute_all_edge_efields
  static void compute_center_efield (Block *block, int dim,
				     std::string center_efield_name,
				     Grouping &prim_group, int stale_depth);

  /// Computes component i of the edge-centered E-field that sits on the faces
  /// of dimensions j and k (i, j, and k are any cyclic permutation of x, y, z).
  /// This uses the component i of the cell-centered E-field and component i
  /// of the E-field on the j and k faces (this comes from fluxes computed for
  /// B-field stored in jflux_group and kflux_group). Additionally, it requires
  /// knowledge of the upwind direction on the j and k faces.
  ///
  /// @param block holds data to be processed
  /// @param dim The component of the edge-centered E-field to compute. Values
  ///     of 0, 1 and 2 correspond to the x, y and z directions, respectively.
  /// @param center_efield_name The name of the field containing component i of
  ///     the cell-centered E-field
  /// @param efield_group This contains the name of the field that will hold
  ///     the values of the computed edge-centered E-field. This must have a
  ///     group called "efield" which contains 3 edge-centered fields (one for
  ///     each component). A given component should be face-centered (should not
  ///     include space for values on the exterior face of the grid) for the
  ///     dimensions that don't match the name of the component and
  ///     cell-centered along the other dimension (e.g. the y-component should
  ///     be face-centered along x and z, but cell-centered along y).
  /// @param jflux_group,kflux_group These contain the fields holding fluxes
  ///     computed along the jth and kth dimensions. The relevant fields are
  ///     face-centered along the j and k dimensions.
  /// @param weight_group holds the temporary weight fields (in a group called
  ///     "weight"). There is a weight field for each spatial direction and
  ///     the field is face-centered along that direction (without including
  ///     values on the exterior faces of the block). The weight fields
  ///     corresponding to j and k indicate the upwind direction along those
  ///     dimensions (this is included to optionally implement the weighting
  ///     scheme used by Athena++ at a later date)
  /// @param stale_depth indicates the current stale_depth for the supplied
  ///     quantities
  ///
  /// @note this function is called in compute_all_edge_efields
  void static compute_edge_efield (Block *block, int dim,
				   std::string center_efield_name,
				   Grouping &efield_group,
				   Grouping &jflux_group,
				   Grouping &kflux_group,
				   Grouping &weight_group,
				   int stale_depth);

  /// Compute the all of the edge-centered electric fields using the current
  /// fluxes and current cell-centered integrable quantities .
  ///
  /// @param block holds data to be processed
  /// @param prim_group holds the fields containing the current values of the
  ///     cell-centered integrable quantities (cell-centered E-field is
  ///     computed from the stored bfield and velocity)
  /// @param xflux_group,yflux_group,zflux_group holds field names where the
  ///     fluxes along the x, y, and z directions are stored. These should
  ///     include fluxes computed for the magnetic fields
  /// @param center_efield_name name of the fields where components of the
  ///     cell-centered electric field can be temporarily stored.
  /// @param efield_group holds field names where the computed edge-centered
  ///     electric fields will be stored. This must have a group called
  ///     "efield" which contains 3 edge-centered fields (one for each
  ///     component). A given component should be face-centered (shouldn't
  ///     include space for values on the exterior faces of the grid) for the
  ///     dimensions that don't match the name of the component and
  ///     cell-centered along the other dimension (e.g. the y-component should
  ///     be face-centered along x and z, but cell-centered along y)
  /// @param weight_group holds the temporary weight fields (in a group called
  ///     "weight"). There is a weight field for each spatial direction and
  ///     the field is face-centered along that direction (without including
  ///     values on the exterior faces of the block). The weight fields for a
  ///     given dimension indicate the upwind direction along that dimension
  ///     (this is included to optionally implement the weighting scheme used
  ///     by Athena++ at a later date)
  /// @param stale_depth the stale depth at the time of this function call
  static void compute_all_edge_efields(Block *block, Grouping &prim_group,
				       Grouping &xflux_group,
				       Grouping &yflux_group,
				       Grouping &zflux_group,
				       std::string center_efield_name,
				       Grouping &efield_group,
				       Grouping &weight_group,
				       int stale_depth);

  /// Updates the face-centered B-field component along the ith dimension using
  /// the jth and kth components of the edge-centered E-field
  ///
  /// @param block holds data to be processed
  /// @param dim The component of the interface B-field to update. Values of 0,
  ///     1 and 2 correspond to the x, y and z directions, respectively.
  /// @param efield_group this contains the edge-centered electric fields used
  ///     to update the interface bfield. This must have a group called
  ///     "efield" which contains 3 edge-centered fields (one for each
  ///     component). A given component should be face-centered (shouldn't
  ///     include space for values on the exterior face of the grid) for the
  ///     dimensions that don't match the name of the component and
  ///     cell-centered along the other dimension (e.g. the y-component
  ///     should be face-centered along x and z, but cell-centered along y).
  /// @param cur_bfieldi_group this contains the current interface B-fields.
  ///     It must have a group called "bfield" holding 3 face-centered fields
  ///     (one for each B-field component). A given component should be
  ///     cell-centered for the dimensions that don't match the name of the
  ///     component and face-centered (including space for values located on
  ///     exterior faces of the block) along the dimension that matches the
  ///     name of the component (e.g. the y-component should be cell-centered
  ///     along x and z, but face-centered along y).
  /// @param out_bfieldi_group this contains the field where the updated
  ///     interface bfield should be stored. The requirements of the contained
  ///     groups are the same as for cur_bfieldi_group. (Note, this can be
  ///     the same Grouping as cur_bfieldi_group.
  /// @param dt The time time-step over which to apply the fluxes
  /// @param stale_depth indicates the current stale_depth for the supplied
  ///     quantities
  static void update_bfield(Block *block, int dim, Grouping &efield_group,
			    Grouping &cur_bfieldi_group,
			    Grouping &out_bfieldi_group,
			    enzo_float dt, int stale_depth);

protected: // attributes

  /// Grouping that holds a group called "bfield". Within that group, there are
  /// three fields. Each field holds a component of the bfield that is
  /// face-centered along the dimension of the component (hence they are
  /// interface bfields)
  Grouping bfieldi_group_;

  /// Contains the relevant current data
  Block* block_;

  /// Grouping of temporary interface bfields identical to all the fields held
  /// by bfieldi_group_. These hold the values of the face-centered fields
  /// computed at the half time-step
  Grouping temp_bfieldi_group_;

  /// Grouping that holds temporary fields inside of a single group called
  /// "weight". "weight" holds three fields, one for each spatial direction.
  /// The weight field for a given direction is face-centered along that
  /// direction (but the field exclude exterior faces of the grid) and keeps
  /// track of the upwind direction.
  Grouping weight_group_;

  /// Grouping of temporary edge-centered fields where the calculated electric
  /// fields get stored. This has one group called "efield" which contains
  /// temporary fields used to store the electric field along the x, y, and z
  /// components. The temporary field will be cell centered along the dimension
  /// of the electric field component and face-centered along the other
  /// components, excluding exterior faces of the mesh (e.g. The x-component
  /// will be cell-centered along the x-direction, but face-centered along the
  /// y- and z- axes)
  Grouping efield_group_;

  /// Holds the name of the temporary field used to store the cell-centered
  /// electric field. This field gets reused for each dimension.
  std::string center_efield_name_;

  /// Number of partial timesteps in a cycle. A value of 1 would means that
  /// there is only a single update over the full timestep. A value of 2 means
  /// that there is a half timestep and a full timestep.
  const int num_partial_timesteps_;

  /// Holds the index of the current partial timestep.
  int partial_timestep_index_;
};
#endif /* ENZO_ENZO_CONSTRAINEDTRANSPORT_HPP */
