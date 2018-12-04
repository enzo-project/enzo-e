
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
				   Grouping &prim_group)
  {
    EnzoBlock * enzo_block = enzo::block(block);
    Field field = enzo_block->data()->field();

    // Load the E-field
    enzo_float* efield = (enzo_float *) field.values(center_efield_id);

    int j = (dim+1)%3;
    int k = (dim+2)%3;

    // Load the jth and kth components of the velocity and cell-centered bfield
    enzo_float* velocity_j = load_grouping_field_(field, &prim_group,
						  "velocity", j);
    enzo_float* velocity_k = load_grouping_field_(field, &prim_group,
						  "velocity", k);
    enzo_float* bfield_j = load_grouping_field_(field, &prim_group,
						"bfield", j);
    enzo_float* bfield_k = load_grouping_field_(field, &prim_group,
						"bfield", k);

    // get interation limits - it includes the ghost zones
    int mx = enzo_block->GridDimension[0];
    int my = enzo_block->GridDimension[1];
    int mz = enzo_block->GridDimension[2];

    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
	for (int ix=0; ix<mx; ix++) {
	  // compute the index
	  int i = ix + mx*(iy + my*iz);
	  efield[i] = velocity_j[i] * bfield_k[i] - velocity_k[i] * bfield_j[i];
	}
      }
    }
  }

  // Computes the edge-centered E-fields pointing in the ith direction
  // It uses the component of the cell-centered E-field pointing in that
  // direction, and the face-centered E-field pointed in that direction
  // the face-centered E-fields are given by elements of jflux_ids and
  // kflux_ids. dim points along i.
  // i, j, and k are any cyclic permutation of x, y, z
  void compute_edge_efield (Block *block, int dim, int efield_id,
			    int center_efield_id, Grouping &jflux_group,
			    Grouping &kflux_group, Grouping &prim_group)
  { }
};

#endif /* ENZO_ENZO_CONSTRAINEDTRANSPORT_HPP */
