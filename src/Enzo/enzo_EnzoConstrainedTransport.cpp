#include "cello.hpp"
#include "enzo.hpp"

void initialize_field_array_(Block *block, EnzoArray<enzo_float> &array,
			     int field_id)
{
  Field field = block->data()->field();
  int mx, my, mz;
  field.dimensions (field_id,&mx,&my,&mz);
  enzo_float *data = (enzo_float *) field.values(field_id);
  array.initialize_wrapper(data, mz, my, mx);
}


void EnzoConstrainedTransport::compute_center_efield (Block *block, int dim,
						      int center_efield_id,
						      Grouping &prim_group)
{
  // Load the E-field
  EnzoArray<enzo_float> efield;
  initialize_field_array_(block, efield, center_efield_id);

  int j = (dim+1)%3;
  int k = (dim+2)%3;

  // Load the jth and kth components of the velocity and cell-centered bfield
  EnzoArray<enzo_float> velocity_j, velocity_k, bfield_j, bfield_k;
  load_grouping_field_(block, prim_group, "velocity", j, velocity_j);
  load_grouping_field_(block, prim_group, "velocity", j, velocity_k);
  load_grouping_field_(block, prim_group, "bfield", j, bfield_j);
  load_grouping_field_(block, prim_group, "bfield", k, bfield_k);

  for (int iz=0; iz<efield.length_dim2(); iz++) {
    for (int iy=0; iy<efield.length_dim1(); iy++) {
      for (int ix=0; ix<efield.length_dim0(); ix++) {
	efield(iz,iy,ix) = (velocity_j(iz,iy,ix) * bfield_k(iz,iy,ix) -
			    velocity_k(iz,iy,ix) * bfield_j(iz,iy,ix));
      }
    }
  }
}

// The following is a helper function that actually computes the component of
// the edge-centered E-field along dimension i.
//
// First, a note on weight arrays:
//   If the i-direction poins along z, we would need to compute
//   dEzdx(ix-1/4,iy-1/2)
//
//   if upwind in positive y direction
//     dEdx(ix-1/4,iy-1/2) = 2.(E(ix,iy-1) - E(ix-1/2,iy-1))/dx
//   if upwind in negative y direction
//     dEdx(ix-1/4,iy-1/2) = 2.(E(ix,iy) - E(ix-1/2,iy))/dx
//   otherwise:
//     dEdx(ix-1/4,iy-1/2) = [(E(ix,iy-1) - E(ix-1/2,iy-1))/dx
//                            + (E(ix,iy) - E(ix-1/2,iy))/dx ]
//   We pulled out a factor of 2 and dx from the derivatives (they cancel)
//  
//   If upwind is in the positive y direction W_y = 1, if downwind W_y = 0,
//     and otherwise W_y = 0.5 (W_y is centered on faces along y-direction)
//
// The form of equations from Stone & Gardiner 09 are as follows. This is
// all necessary to compute Ez(ix-1/2, iy-1/2, iz): 
//  dEdj_h = dEdx(ix-1/4,iy-1/2) =
//       W_y(    ix,iy-1/2)  * (E(    ix,  iy-1) - E(ix-1/2,  iy-1)) +
//    (1-W_y(    ix,iy-1/2)) * (E(    ix,    iy) - E(ix-1/2,    iy))
//  dEdj_l = dEdx(ix-3/4,iy-1/2) =
//       W_y(  ix-1,iy-1/2)  * (E(ix-1/2,  iy-1) - E(  ix-1,  iy-1)) +
//    (1-W_y(  ix-1,iy-1/2)) * (E(ix-1/2,    iy) - E(  ix-1,    iy))
//  dEdk_h = dEdy(ix-1/2,iy-1/4) =
//       W_x(ix-1/2,   iy)  * (E(   ix-1,    iy) - E(  ix-1,iy-1/2)) +
//    (1-W_x(ix-1/2,   iy)) * (E(     ix,    iy) - E(    ix,iy-1/2))
//  dEdk_l = dEdy(ix-1/2,iy-3/4) =
//       W_x(ix-1/2, iy-1)  * (E(   ix-1,iy-1/2) - E(  ix-1,  iy-1)) +
//    (1-W_x(ix-1/2, iy-1)) * (E(     ix,iy-1/2) - E(    ix,  iy-1))
//
// This solution to the problem is easier to program if we reframe it as the
// calculation of Ez(ix+1/2, iy+1/2, iz). The above equations are then
// rewritten as:
//  dEdj_h = dEdx(ix+3/4,iy+1/2) =
//       W_y(  ix+1,iy+1/2)  * (E(  ix+1,    iy) - E(ix+1/2,    iy)) +
//    (1-W_y(  ix+1,iy+1/2)) * (E(  ix+1,  iy+1) - E(ix+1/2,  iy+1))
//  dEdj_l = dEdx(ix+1/4,iy+1/2) =
//       W_y(    ix,iy+1/2)  * (E(ix+1/2,    iy) - E(    ix,    iy)) +
//    (1-W_y(    ix,iy+1/2)) * (E(ix+1/2,  iy+1) - E(    ix,  iy+1))
//  dEdk_h = dEdy(ix+1/2,iy+3/4) =
//       W_x(ix+1/2,  iy+1)  * (E(     ix,  iy+1) - E(    ix,iy+1/2)) +
//    (1-W_x(ix+1/2,  iy+1)) * (E(   ix+1,  iy+1) - E(  ix+1,iy+1/2))
//  dEdk_l = dEdy(ix+1/2,iy+1/4) =
//       W_x(ix+1/2,   iy)  * (E(     ix,iy+1/2) - E(    ix,    iy)) +
//    (1-W_x(ix+1/2,   iy)) * (E(   ix+1,iy+1/2) - E(  ix+1,    iy))
//
// Now let's rewrite in totally Generalized notation:  [z->i, x->j, y->k] 
//  dEdj_h = dEdj(k+1/2,j+3/4,i) =
//       W_k(k+1/2,  j+1)  * (E(    k,  j+1) - E(    k,j+1/2)) +
//    (1-W_k(k+1/2,  j+1)) * (E(  k+1,  j+1) - E(  k+1,j+1/2))
//  dEdj_l = dEdj(k+1/2,j+1/4,i) =
//       W_k(k+1/2,    j)  * (E(    k,j+1/2) - E(    k,    j)) +
//    (1-W_k(k+1/2,    j)) * (E(  k+1,j+1/2) - E(  k+1,    j))
//  dEdk_h = dEdk(k+3/4,j+1/2,i) =
//       W_j(  k+1,j+1/2)  * (E(  k+1,    j) - E(k+1/2,    j)) +
//    (1-W_j(  k+1,j+1/2)) * (E(  k+1,  j+1) - E(k+1/2,  j+1))
//  dEdk_l = dEdk(k+1/4,j+1/2,i) =
//       W_j(    k,j+1/2)  * (E(k+1/2,     j) - E(    k,    j)) +
//    (1-W_j(    k,j+1/2)) * (E(k+1/2,   j+1) - E(    k,  j+1))
//
//
// Define subarrays of W_j, W_k, E_cen, E_j, and E_k
//   Weights:
//     Wj(k,j,i)       =   W_j(    k,j+1/2,i)
//     Wj_kp1(k,j,i)   =   W_j(  k+1,j+1/2,i)
//     Wk(k,j,i)       =   W_k(k+1/2,    j,i)
//     Wk_jp1(k,j,i)   =   W_k(k+1/2,  j+1,i)
//
//   Cell-Centered Efield (i-component):
//     Ec(k,j,i)       =     E(    k,    j,i)
//     Ec_jp1(k,j,i)   =     E(    1,  j+1,i)
//     Ec_kp1(k,j,i)   =     E(  k+1,    j,i)
//     Ec_jkp1(k,j,i)  =     E(  k+1,  j+1,i)
//
//   face-centered E-field (i-component):
//     Ej(k,j,i)       =     E(    k,j+1/2,i)
//     Ej_kp1(k,j,i)   =     E(  k+1,j+1/2,i)
//     Ek(k,j,i)       =     E(k+1/2,    j,i)
//     Ek_jp1(k,j,i)   =     E(k+1/2,  j+1,i)
//
//   edge-centered E-field (i-component):
//     Eedge(k,j,i)    =     E(k+1/2,j+1/2,i)
// To be independent of reconstruction method, compute E-field at all edges
// that lie within the mesh.

void compute_edge_(int xstart, int ystart, int zstart,
		   int xstop, int ystop, int zstop,
		   EnzoArray<enzo_float> &Eedge,
		   EnzoArray<enzo_float> &Wj, EnzoArray<enzo_float> &Wj_kp1,
		   EnzoArray<enzo_float> &Wk, EnzoArray<enzo_float> &Wk_jp1,
		   EnzoArray<enzo_float> &Ec, EnzoArray<enzo_float> &Ec_jkp1,
		   EnzoArray<enzo_float> &Ec_jp1, EnzoArray<enzo_float> &Ec_kp1,
		   EnzoArray<enzo_float> &Ej, EnzoArray<enzo_float> &Ej_kp1,
		   EnzoArray<enzo_float> &Ek, EnzoArray<enzo_float> &Ek_jp1)
{
  for (int iz = zstart; iz < zstop; iz++){
    for (int iy = ystart; iy < ystop; iy++){
      for (int ix = xstart; ix < xstop; ix++){

	enzo_float dEdj_l, dEdj_h, dEdk_l, dEdk_h;

	//  dEdj(k+1/2,j+3/4,i) =
	//       W_k(k+1/2,  j+1)  * (E(    k,  j+1) - E(    k,j+1/2)) +
	//    (1-W_k(k+1/2,  j+1)) * (E(  k+1,  j+1) - E(  k+1,j+1/2))
	dEdj_h =
	       Wk_jp1(iz,iy,ix)  * ( Ec_jp1(iz,iy,ix) -     Ej(iz,iy,ix)) +
	  (1 - Wk_jp1(iz,iy,ix)) * (Ec_jkp1(iz,iy,ix) - Ej_kp1(iz,iy,ix));

	//  dEdj(k+1/2,j+1/4,i) =
	//       W_k(k+1/2,    j)  * (E(    k,j+1/2) - E(    k,    j)) +
	//    (1-W_k(k+1/2,    j)) * (E(  k+1,j+1/2) - E(  k+1,    j))
	dEdj_l =
	           Wk(iz,iy,ix)  * (     Ej(iz,iy,ix) -     Ec(iz,iy,ix)) +
	  (1 -     Wk(iz,iy,ix)) * ( Ej_kp1(iz,iy,ix) - Ec_kp1(iz,iy,ix));

	//  dEdk(k+3/4,j+1/2,i) =
	//       W_j(  k+1,j+1/2)  * (E(  k+1,    j) - E(k+1/2,    j)) +
	//    (1-W_j(  k+1,j+1/2)) * (E(  k+1,  j+1) - E(k+1/2,  j+1))
	dEdk_h =
	       Wj_kp1(iz,iy,ix)  * ( Ec_kp1(iz,iy,ix) -     Ek(iz,iy,ix)) +
	  (1 - Wj_kp1(iz,iy,ix)) * (Ec_jkp1(iz,iy,ix) - Ek_jp1(iz,iy,ix));

	//  dEdk(k+1/4,j+1/2,i) =
	//       W_j(    k,j+1/2)  * (E(k+1/2,     j) - E(    k,    j)) +
	//    (1-W_j(    k,j+1/2)) * (E(k+1/2,   j+1) - E(    k,  j+1))
	dEdk_l =
	           Wj(iz,iy,ix)  * (     Ek(iz,iy,ix) -     Ec(iz,iy,ix)) +
	  (1 -     Wj(iz,iy,ix)) * ( Ek_jp1(iz,iy,ix) - Ec_jp1(iz,iy,ix));

	Eedge(iz,iy,ix) = 0.25*(Ej(iz,iy,ix) + Ej_kp1(iz,iy,ix) +
				Ek(iz,iy,ix) + Ek_jp1(iz,iy,ix) +
				(dEdj_h - dEdj_l) + (dEdk_h - dEdk_l));		
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
//
// This Method is applicable for:
//    - 2D array (z-axis only has 1 entry), with dim = 2
// 
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
						    int center_efield_id,
						    Grouping &efield_group,
						    Grouping &jflux_group,
						    Grouping &kflux_group,
						    Grouping &prim_group,
						    Grouping &weight_group)
{
  // determine alignment of j,k axes with respect to x,y, and z
  int dxdj = 0; int dydj = 0; int dzdj=0;
  int dxdk = 0; int dydk = 0; int dzdk=0;
  if (dim == 0){
    dydj = 1;
    dzdk = 1;
  } else if (dim == 1){
    dzdj = 1;
    dxdk = 1;
  } else {
    dxdj = 1;
    dydk = 1;
  }

  // Initialize Cell-Centered E-fields
  EnzoArray<enzo_float> Ec, Ec_jp1, Ec_kp1, Ec_jkp1;
  initialize_field_array_(block, Ec, center_efield_id);
  Ec_jp1.initialize_subarray(Ec, dzdj, Ec.length_dim2(),  //zstart : zstop
			         dydj, Ec.length_dim1(),  //ystart : ystop
			         dxdj, Ec.length_dim0()); //xstart : xstop
  Ec_kp1.initialize_subarray(Ec, dzdk, Ec.length_dim2(),
			         dydk, Ec.length_dim1(),
			         dxdk, Ec.length_dim0());
  Ec_jkp1.initialize_subarray(Ec, dzdj+dzdk, Ec.length_dim2(),
			          dydj+dydk, Ec.length_dim1(),
			          dxdj+dxdk, Ec.length_dim0());

  // Initialize edge-centered Efield [it maps (k,j,i) -> (k+1/2,j+1/2,i)]
  EnzoArray<enzo_float> Eedge;
  load_temp_interface_grouping_field_(block, efield_group, "efield", dim,
				      Eedge, dim == 0, dim == 1, dim == 2);

  // Initialize face-centered E-fields
  EnzoArray<enzo_float> Ej, Ej_kp1, Ek, Ek_jp1;
  load_temp_interface_grouping_field_(block, jflux_group, "bfield", dim,
				      Ej, dxdj==0, dydj==0, dzdj==0);
  Ej_kp1.initialize_subarray(Ej, dzdk, Ej.length_dim2(),
			         dydk, Ej.length_dim1(),
			         dxdk, Ej.length_dim0());
  load_temp_interface_grouping_field_(block, kflux_group, "bfield", dim,
				      Ek, dxdk==0, dydk==0, dzdk==0);
  Ek_jp1.initialize_subarray(Ek, dzdj, Ek.length_dim2(),
			         dydj, Ek.length_dim1(),
			         dxdj, Ek.length_dim0());

  // Initialize the weight arrays
  EnzoArray<enzo_float> Wj, Wj_kp1, Wk, Wk_jp1;
  load_temp_interface_grouping_field_(block, weight_group, "weight", (dim+1)%3,
				      Wj, dxdj==0, dydj==0, dzdj==0);
  Wj_kp1.initialize_subarray(Wj, dzdk, Wj.length_dim2(),
			         dydk, Wj.length_dim1(),
			         dxdk, Wj.length_dim0());
  load_temp_interface_grouping_field_(block, weight_group, "weight", (dim+2)%3,
				      Wk, dxdk==0, dydk==0, dzdk==0);
  Wk_jp1.initialize_subarray(Wj, dzdj, Wk.length_dim2(),
			         dydj, Wk.length_dim1(),
			         dxdj, Wk.length_dim0());

  // Integration limits
  //
  // If computing the edge E-field along z-direction:
  //    - If grid is 2D (the grid only has 1 cell along z-component),
  //      need to iterate from 0 to 1 for z. (Will never compute E-field
  //      along any other dimension for a 2D grid).
  //    - If grid is 3D (the grid has more than 1 cell along the z-component),
  //      there is no need to compute the e-field at iz=0 or iz= mz-1 since
  //      we have E-fields on the exterior of the mesh. (This logic applies to
  //      components of other dimensions). Generalizing to computing
  //      E-field along dimension i, then we need to compute it at i=1 up to
  //      (but not including imax-1)
  // For all cases, if we are computing the E-field along dimension i, then we
  // need to compute it at:  j = 1/2 up to (but not including) j = jmax-1/2
  //                         k = 1/2 up to (but not including) k = kmax-1/2
  //
  // Note if an quantitiy is face-centered along dimension dim:
  //    idim maps to idim+1/2
  //
  // To summarize:
  //    2D grid with i aligned along z:  istart = 0     istop = imax = 1
  //    3D grid:                         istart = 1     istop = imax - 1
  //    In all cases:                    jstart = 0     jstop = jmax - 1
  //                                     kstart = 0     kstop = kmax - 1

  int xstart = 1 - dxdj - dxdk; // if dim==0: 1, otherwise: 0
  int ystart = 1 - dydj - dydk;
  int zstart = (Ec.length_dim2() == 1) ? 0 : 1 - dzdj - dzdk;

  int xstop = Ec.length_dim0() - 1;
  int ystop = Ec.length_dim1() - 1;
  int zstop = (Ec.length_dim2() == 1) ? 1 : Ec.length_dim2() - 1;

  compute_edge_(xstart, ystart, zstart, xstop, ystop, zstop,
		Eedge, Wj, Wj_kp1, Wk, Wk_jp1, Ec, Ec_jkp1, Ec_jp1, Ec_kp1,
		Ej, Ej_kp1, Ek, Ek_jp1);

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
