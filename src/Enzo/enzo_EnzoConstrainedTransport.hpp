
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

  void compute_center_efield (Block *block, int dim, int center_efield_id,
			      Grouping &prim_group);
  
  // Computes the edge-centered E-fields pointing in the ith direction
  // It uses the component of the cell-centered E-field pointing in that
  // direction, and the face-centered E-field pointed in that direction
  // the face-centered E-fields are given by elements of jflux_ids and
  // kflux_ids. dim points along i.
  // i, j, and k are any cyclic permutation of x, y, z
  void compute_edge_efield (Block *block, int dim, int center_efield_id,
			    Grouping &efield_group, Grouping &jflux_group,
			    Grouping &kflux_group, Grouping &prim_group,
			    Grouping &weight_group);

  void update_bfield(Block *block, int dim, Grouping &efield_group,
		     Grouping &cur_bfieldi_group, Grouping &out_bfieldi_group,
		     enzo_float dt);

  // if compute_outer is false, then the centered B-field is computed for all
  // cells in the grid except at i = 0 and i = imax-1, where i is the index for
  // the axis aligned with dim. If true, then computes the values at all cells,
  // including the cells at i = 0 and i = imax-1. If ture, bfieldi_group is
  // assumed to include interface bfields with values for the exterior of the
  // grid.
  //
  // compute_outer should only be set to true when calculating the cell-
  // centered magnetic fields for the first time. During the updates from
  // Constrained Transport, it should always be false ( since the temporary
  // bfields at the half time-step do not include values for the exterior
  // faces of the grid).
  void compute_center_bfield(Block *block, int dim, Grouping &cons_group,
			     Grouping &bfieldi_group, bool compute_outer);

};
#endif /* ENZO_ENZO_CONSTRAINEDTRANSPORT_HPP */
