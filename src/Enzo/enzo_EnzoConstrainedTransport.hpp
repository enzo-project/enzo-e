
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

  void compute_cell_center_efield (Block *block, int dim, int center_efield_id,
				   std::vector<int> &prim_ids)
  { }

  // Computes the edge-centered E-fields pointing in the ith direction
  // It uses the component of the cell-centered E-field pointing in that
  // direction, and the face-centered E-field pointed in that direction
  // the face-centered E-fields are given by elements of jflux_ids and
  // kflux_ids. dim points along i.
  // i, j, and k are any cyclic permutation of x, y, z
  void compute_edge_efield (Block *block, int dim, int efield_id,
			    int center_efield_id,
			    std::vector<int> &jflux_ids,
			    std::vector<int> &kflux_ids,
			    std::vector<int> &prim_ids)
  { }
};

#endif /* ENZO_ENZO_CONSTRAINEDTRANSPORT_HPP */
