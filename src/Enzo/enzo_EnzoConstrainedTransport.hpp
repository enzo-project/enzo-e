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

public: // interface

  /// Create a new EnzoConstrainedTransport
  EnzoConstrainedTransport() throw()
  {}

  /// Virtual destructor
  virtual ~EnzoConstrainedTransport()
  {  }

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
  void compute_center_efield (Block *block, int dim,
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
  void compute_edge_efield (Block *block, int dim,
			    std::string center_efield_name,
			    Grouping &efield_group, Grouping &jflux_group,
			    Grouping &kflux_group, Grouping &weight_group,
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
  ///
  void compute_all_edge_efields(Block *block, Grouping &prim_group,
				Grouping &xflux_group, Grouping &yflux_group,
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
  void update_bfield(Block *block, int dim, Grouping &efield_group,
		     Grouping &cur_bfieldi_group, Grouping &out_bfieldi_group,
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
  void compute_center_bfield(Block *block, int dim, Grouping &bfieldc_group,
			     Grouping &bfieldi_group, int stale_depth = 0);

};
#endif /* ENZO_ENZO_CONSTRAINEDTRANSPORT_HPP */
