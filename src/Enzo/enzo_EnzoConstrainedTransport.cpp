#include "cello.hpp"
#include "enzo.hpp"

void EnzoConstrainedTransport::compute_center_efield (Block *block, int dim,
						      int center_efield_id,
						      Grouping &prim_group)
{
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // Load the E-field
  enzo_float* efield = (enzo_float *) field.values(center_efield_id);

  int j = (dim+1)%3;
  int k = (dim+2)%3;

  // Load the jth and kth components of the velocity and cell-centered bfield
  enzo_float* velocity_j = load_grouping_field_(&field, &prim_group,
						"velocity", j);
  enzo_float* velocity_k = load_grouping_field_(&field, &prim_group,
						"velocity", k);
  enzo_float* bfield_j = load_grouping_field_(&field, &prim_group,
					      "bfield", j);
  enzo_float* bfield_k = load_grouping_field_(&field, &prim_group,
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


// helper function that computes sizes and offsets of the arrays
// used to compute the edge-centered E-fields
void compute_sizes_and_offsets_(int dim, int mx, int my, int mz,
				int &ec_mz, int &ec_my, int &ec_mx,
				int &jfc_mz, int &jfc_my, int &jfc_mx,
				int &kfc_mz, int &kfc_my, int &kfc_mx,
				int &c_off_j, int &c_off_k,
				int &jflux_off_k, int &kflux_off_j)
{
  // dim is an integer corresponding to the ith direction
  // -- ec_mz, ec_my, ec_mx are the dimensions of an array of values
  //    centered on edges. Specifically, the values are on the faces of
  //    dimensions j and k
  // -- jfc_mz, jfc_my, jfc_mx are the dimensions of an array of values
  //    centered on the faces of the mesh along dimension j
  // -- kfc_mz, kfc_my, kfc_mx are the dimensions of an array of values
  //    centered on the faces of the mesh along dimension k
  // -- c_off_j and c_off_k are the offsets for the cell centered field to
  //    access the next value in the positive j and k directions.
  // -- jflux_off_k is the offset (in array index) required to access
  //    the next value of the flux field, centered on faces along dimension j,
  //    in the positive k direction.
  // -- kflux_off_j is the offset (in array index) required to access
  //    the next value of the flux field, centered on faces along dimension k,
  //    in the positive j direction.

  ec_mx = mx; ec_my = my; ec_mz = mz;
  jfc_mx = mx; jfc_my = my; jfc_mz = mz;
  kfc_mx = mx; kfc_my = my; kfc_mz = mz;
  
  if (dim == 0) {
    ec_my++;
    ec_mz++;

    jfc_my++;
    kfc_mz++;

    c_off_j = mx;
    c_off_k = mx*my;

    jflux_off_k = jfc_mx*jfc_my;
    kflux_off_j = kfc_mx;    
  } else if (dim == 1) {
    ec_mx++;
    ec_mz++;

    jfc_mz++;
    kfc_mz++;

    c_off_j = mx*my;
    c_off_k = 1;

    jflux_off_k = 1;
    kflux_off_j = kfc_mx*kfc_my;
  } else {
    ec_mx++;
    ec_my++;

    jfc_mx++;
    kfc_my++;

    c_off_j = 1;
    c_off_k = mx;

    jflux_off_k = jfc_mx;
    kflux_off_j = 1;
  }
}


enzo_float compute_edge_(int c_ind, int jfc_ind, int kfc_ind, int c_off_j,
			 int c_off_k, int jfc_off_k, int kfc_off_j,
			 enzo_float* E_cen, enzo_float* E_j,
			 enzo_float* E_k, enzo_float* W_j, enzo_float* W_k)
{
  // helper function that actually computes the component of the edge-centered
  // E-field along dimension i at location (i,j-1/2,k-1/2)
  
  // All E-field components used in this function point in the i-direction
  // -- c_ind is the index of cell-centered fields that corresponds to location
  //    (i, j, k)
  // -- jfc_ind/kfc_ind is the index of the fields centered in on the faces of
  //    dimension j/k that corresponds to location (i, j-1/2,k)/(i,j,k-1/2)
  // -- c_off_j/c_off_k is the offset in the array index required to incrament
  //    the j/k index of the cell-centered field
  // -- jfc_off_k is the offset (in the array index) required to access
  //    the next value of the field, centered on faces along dimension j,
  //    in the positive k direction (e.g. jfc_ind-jfc_off_k corresponds to
  //    the location (i, j-1/2,k-1))
  // -- kfc_off_j is the offset (in the array index) required to access
  //    the next value of the field, centered on faces along dimension k,
  //    in the positive j direction (e.g. kfc_ind-kfc_off_j corresponds to
  //    the location (i, j-1,k-1/2))
  // -- E_center corresponds to the cell-centered E-field component
  // -- E_j/E_k corresponds to the E-field component centered on faces along
  //    the j/k direction
  // -- w_j/w_k corresponds to the weight field centered on faces along the j/k
  //    direction. It can have values of 1, 0.5, and 0 (indicating that it can
  //    the upwind direction is to the positive j/k direction, there is no
  //    upwind direction, or the upwind direction is to the negative j/k
  //    direction). These weight values can use the same more sophisticated
  //    weighting scheme as Athena++

  enzo_float dEdj_l, dEdj_h, dEdk_l, dEdk_h;

  // If the i-direction poins along z:
  //  Computing dEzdx(ix-1/4,iy-1/2)
  //  if upwind in positive y direction
  //    dEdx(ix-1/4,iy-1/2) = 2.(E(ix,iy-1) - E(ix-1/2,iy-1))/dx
  //  if upwind in negative y direction
  //    dEdx(ix-1/4,iy-1/2) = 2.(E(ix,iy) - E(ix-1/2,iy))/dx
  //  otherwise:
  //    dEdx(ix-1/4,iy-1/2) = [(E(ix,iy-1) - E(ix-1/2,iy-1))/dx
  //                           + (E(ix,iy) - E(ix-1/2,iy))/dx ]
  //  We pull out a factor of 2 and dx from the derivatives (they cancel)
  //
  //  If upwind is in the positive y direction W_y = 1, if downwind W_y = 0,
  //    and otherwise W_y = 0.5 (W_y is centered on faces along y-direction)

  //  dEdj_h = dEdx(ix-1/4,iy-1/2) =
  //       W_y(    ix,iy-1/2)  * (E(    ix,  iy-1) - E(ix-1/2,  iy-1)) +
  //    (1-W_y(    ix,iy-1/2)) * (E(    ix,    iy) - E(ix-1/2,    iy))
  dEdj_h = W_k[kfc_ind]*(E_cen[c_ind-c_off_k]                      //j,k-1
			 - E_j[jfc_ind-jfc_off_k])                 //j-1/2,k-1
    + (1.-W_k[kfc_ind])*(E_cen[c_ind]                              //j,k
			 - E_j[jfc_ind]);                          //j-1/2,k

  //  dEdj_l = dEdx(ix-3/4,iy-1/2) =
  //       W_y(  ix-1,iy-1/2)  * (E(ix-1/2,  iy-1) - E(  ix-1,  iy-1)) +
  //    (1-W_y(  ix-1,iy-1/2)) * (E(ix-1/2,    iy) - E(  ix-1,    iy))
  dEdj_l = W_k[kfc_ind-kfc_off_j]*(E_j[jfc_ind-jfc_off_k]          //j-1/2,k-1
				   - E_cen[c_ind-c_off_j-c_off_k]) //j-1,k-1
    + (1.-W_k[kfc_ind-kfc_off_j])*(E_j[jfc_ind]                    //j-1/2,k
				  - E_cen[c_ind-c_off_j]);         //j-1,k

  //  dEdk_h = dEdy(ix-1/2,iy-1/4) =
  //       W_x(ix-1/2,   iy)  * (E(   ix-1,    iy) - E(  ix-1,iy-1/2)) +
  //    (1-W_x(ix-1/2,   iy)) * (E(     ix,    iy) - E(    ix,iy-1/2))
  dEdk_h = W_j[jfc_ind]*(E_cen[c_ind-c_off_j]                      //j-1,k
			 - E_k[kfc_ind-kfc_off_j])                 //j-1,k-1/2
    + (1.-W_j[jfc_ind])*(E_cen[c_ind]                              //j,k
			 - E_k[kfc_ind]);                          //j,k-1/2

  //  dEdk_l = dEdy(ix-1/2,iy-3/4) =
  //       W_x(ix-1/2, iy-1)  * (E(   ix-1,iy-1/2) - E(  ix-1,  iy-1)) +
  //    (1-W_x(ix-1/2, iy-1)) * (E(     ix,iy-1/2) - E(    ix,  iy-1))
  dEdk_l = W_j[jfc_ind-jfc_off_k]*(E_k[kfc_ind-kfc_off_j]          //j-1,k-1/2
				   - E_cen[c_ind-c_off_j-c_off_k]) //j-1,k-1
    + (1.-W_j[jfc_ind-jfc_off_k])*(E_k[kfc_ind]                    //j,k-1/2
				  - E_cen[c_ind-c_off_k]);         //j,k-1

  return 0.25*((E_j[jfc_ind] + E_j[jfc_ind-jfc_off_k] +
		E_k[kfc_ind] + E_k[kfc_ind-kfc_off_j]) +
	       (dEdj_h - dEdj_l) + (dEdk_h - dEdk_l));
}

// Computes the edge-centered E-fields pointing in the ith direction
// It uses the component of the cell-centered E-field pointing in that
// direction, and the face-centered E-field pointed in that direction
// the face-centered E-fields are given by elements of jflux_ids and
// kflux_ids. dim points along i.
// i, j, and k are any cyclic permutation of x, y, z
//
// Not thrilled with how complicated this is, may rewrite the function more
// explicitly (having code devoted to different directions)
//
// The Athena++ code calculates a quantity they refer to as v_over_c at all
// cell-faces.
//   - Basically, this tells them how much to weight derivatives while
//     computing the edge_efield. If deemed necessary, these values can be
//     passed as weight_group
//   - Weight_group includes 3 temporary fields centered on the faces looking
//     down the x, y, and z direction. Currently, it expects values of 1 and 0
//     to indicate that the upwind direction is in positive and negative
//     direction, or 0 to indicate no upwind direction.
void EnzoConstrainedTransport::compute_edge_efield (Block *block, int dim,
						    int efield_id,
						    int center_efield_id,
						    Grouping &jflux_group,
						    Grouping &kflux_group,
						    Grouping &prim_group,
						    Grouping &weight_group)
{
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // Load the E-fields - only looking at the component along i-direction
  enzo_float* E_edge = (enzo_float *) field.values(efield_id);
  // E_edge points along i dimension with shape (nk+1, nj+1, ni)
  // technically, only need to track (nk-1,nj-1,ni)
  enzo_float* E_center = (enzo_float *) field.values(center_efield_id);
  enzo_float *E_jface = load_grouping_field_(&field, &jflux_group, "bfield",
					     dim);
  enzo_float *E_kface = load_grouping_field_(&field, &kflux_group, "bfield",
					     dim);
  // Load the weights
  enzo_float *W_jface = load_grouping_field_(&field, &weight_group, "weight",
					     (dim+1)%3);
  enzo_float *W_kface = load_grouping_field_(&field, &weight_group, "weight",
					     (dim+2)%3);

  // cell-centered iteration dimensions
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];

  // compute edge-centered dimensions and offset between neigboring elements
  // along j and k dimensions for face-centered and cell-centered quantites
  // (the initial values at declaration don't really matter)
  int ec_mx = mx; int ec_my = my; int ec_mz = mz;
  int jfc_mx = mx; int jfc_my = my; int jfc_mz = mz;
  int kfc_mx = mx; int kfc_my = my; int kfc_mz = mz;
  int c_off_j=0; int c_off_k = 0; int jfc_off_k =0; int kfc_off_j = 0;

  compute_sizes_and_offsets_(dim, mx, my, mz, ec_mz, ec_my, ec_mx,
			     jfc_mz, jfc_my, jfc_mx,
			     kfc_mz, kfc_my, kfc_mx,
			     c_off_j, c_off_k,
			     jfc_off_k, kfc_off_j);

  // If computing the edge E-field along z-direction, there is no point in
  // computing it anywhere with iz=1/2 or iz=mz-1/2 since it wouldn't be used
  // for anything. The relevant centered B-fields would also requre
  // edge-centered E-fields on the exterior of the mesh. The same logic applies
  // for other dimensions
  // To be independent of reconstruction method, compute E-field at all edges
  // that lie within the mesh.

  for (int iz=1; iz<mz; iz++) {
    for (int iy=1; iy<my; iy++) {
      for (int ix=1; ix<mx; ix++) {
	  
	// for edge-centered quantities
	//  -dim = 2: iz,iy,ix corresponds to (    iz, iy-1/2, ix-1/2)
	//  -dim = 1: iz,iy,ix corresponds to (iz-1/2,     iy, ix-1/2)
	//  -dim = 0: iz,iy,ix corresponds to (iz-1/2, iy-1/2,     ix)

	// face-centered quantities - centered along
	//  -dim = 2: iz,iy,ix corresponds to (iz-1/2,     iy,     ix)
	//  -dim = 1: iz,iy,ix corresponds to (    iz, iy-1/2,     ix)
	//  -dim = 0: iz,iy,ix corresponds to (    iz,     iy, ix-1/2)

	// compute the indices of the face-centered, cell-centered, and
	// edge-centered fields
	int jfc_ind = ix + jfc_mx*(iy + jfc_my*iz);
	int kfc_ind = ix + kfc_mx*(iy + kfc_my*iz);
	int c_ind = ix + mx*(iy + my*iz);
	int ec_ind = ix + ec_mx*(iy + ec_my*iz);

	E_edge[ec_ind] = compute_edge_(c_ind, jfc_ind, kfc_ind, c_off_j,
				       c_off_k, jfc_off_k, kfc_off_j, E_center,
				       E_jface, E_kface, W_jface, W_kface);
      }
    }
  }
};

void EnzoConstrainedTransport::update_bfield(Block *block, int dim,
					     Grouping &efield_group,
					     Grouping &cur_cons_group,
					     Grouping &out_cons_group,
					     enzo_float dt)
{
  // computing the face-centered B-field component along the ith dimension
  Field field = block->data()->field();
  EnzoBlock * enzo_block = enzo::block(block);

  // cell-centered iteration dimensions
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];

  // load the fields
  enzo_float *cur_bfield = load_grouping_field_(&field, &cur_cons_group,
						"bfieldi", dim);
  enzo_float *out_bfield = load_grouping_field_(&field, &out_cons_group,
						"bfieldi", dim);
  // Load e-fields
  enzo_float *E_j = load_grouping_field_(&field, &efield_group,
					 "efield", (dim+1)%3);
  enzo_float *E_k = load_grouping_field_(&field, &efield_group,
					 "efield", (dim+2)%3);

  // widths of cells along j and k directions
  enzo_float dtdj = dt/enzo_block->CellWidth[(dim+1)%3];
  enzo_float dtdk = dt/enzo_block->CellWidth[(dim+2)%3];

  // face-centered and edge-centered dimensions and edge-centered offsets
  int fc_mx = mx; int fc_my = my; int fc_mz = mz;
  // all edge-centered fields share one edge with the face in i direction
  int ecij_mx = mx+1; int ecij_my = my + 1; int ecij_mz = mz + 1;
  int ecik_mx = mx+1; int ecik_my = my + 1; int ecik_mz = mz + 1;
  // the offset of the edge-centered efield on the edge of the j/k direction
  // shifts the j/k index
  int ecij_offset, ecik_offset;

  if (dim == 0){
    // dim points in x-direction
    fc_mx++; ecij_mz--; ecik_my--;
    ecij_offset = ecij_mx;
    ecik_offset = ecik_mx * ecik_my;
  } else if (dim == 1){
    fc_my++; ecij_mx--; ecik_mz--;
    ecij_offset = ecij_mx*ecij_my;
    ecik_offset = 1;
  } else {
    fc_mz++; ecij_my--; ecik_mx--;
    ecij_offset = 1;
    ecik_offset = ecik_mx;
  }

  for (int iz=1; iz<fc_mz-1; iz++) {
    for (int iy=1; iy<fc_my-1; iy++) {
      for (int ix=1; ix<fc_mx-1; ix++) {
	// bfield is centered on face in the i direction
	int fc_ind = ix + fc_mx*(iy + fc_my*iz);

	// E-field in k direction centered on edges joining the i and j faces
	int ecij_ind = ix + ecij_mx*(iy + ecij_my*iz);

	// E-field in j direction centered on edges joining the i and k faces
	int ecik_ind = ix + ecik_mx*(iy + ecik_my*iz);

	// if dim were equal to zero, then we would compute:
	//   Bnew_x(i-1/2, j, k) = Bold_x(i-1/2,j,k)
	//     - dt/dy*(E_z(i-1/2,j+1/2,    k) - E_z(i-1/2,j-1/2,    k))
	//     + dt/dz*(E_y(i-1/2,    j,k+1/2) - E_z(i-1/2,    j,k-1/2))

	out_bfield[fc_ind] = cur_bfield[fc_ind]
	  - dtdj*(E_k[ecij_ind + ecij_offset] - E_k[ecij_ind])
	  + dtdk*(E_j[ecik_ind + ecik_offset] - E_j[ecik_ind]);
      }
    }
  }
}

void EnzoConstrainedTransport::compute_center_bfield(Block *block, int dim,
						     Grouping &cons_group,
						     enzo_float dt)
{
  Field field = block->data()->field();

  EnzoBlock * enzo_block = enzo::block(block);

  // cell-centered iteration dimensions
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];

  // compute face-centered dimensions and offset between neigboring elements
  // along dim (offset is the same for cell-centered and face-centered arrays)
  int fc_mx = mx; int fc_my = my; int fc_mz = mz;

  int offset;
  if (dim == 0) {
    fc_mx++;
    offset = 1;
  } else if (dim == 1) {
    fc_my++;
    offset = mx;
  } else {
    fc_mz++;
    offset = mx*my;
  }

  enzo_float *bfield = load_grouping_field_(&field, &cons_group,
					    "bfieldi", dim);
  enzo_float *bfield_center = load_grouping_field_(&field, &cons_group,
						   "bfield", dim);

  for (int iz=1; iz<mz-1; iz++) {
    for (int iy=1; iy<my-1; iy++) {
      for (int ix=1; ix<mx-1; ix++) {
	// compute the index of the cell-centered and face-centered
	int fc_ind = ix + fc_mx*(iy + fc_my*iz);    
	int c_ind = ix + mx*(iy + my*iz);

	bfield_center[c_ind] = 0.5*(bfield[fc_ind+offset] + bfield[fc_ind]);
      }
    }
  }
}
