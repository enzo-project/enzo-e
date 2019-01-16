
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

  void compute_center_bfield(Block *block, int dim, Grouping &cons_group,
			     Grouping &bfieldi_group);

};
#endif /* ENZO_ENZO_CONSTRAINEDTRANSPORT_HPP */
