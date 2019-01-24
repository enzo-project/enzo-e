#include "cello.hpp"
#include "enzo.hpp"

void EnzoConstrainedTransport::compute_center_efield (Block *block, int dim,
						      int center_efield_id,
						      Grouping &prim_group)
{
  // Load the E-field
  EnzoArray<enzo_float> efield;
  EnzoFieldArrayFactory array_factory;
  array_factory.initialize_field_array(block, efield, center_efield_id);

  int j = (dim+1)%3;
  int k = (dim+2)%3;

  // Load the jth and kth components of the velocity and cell-centered bfield
  EnzoArray<enzo_float> velocity_j, velocity_k, bfield_j, bfield_k;
  array_factory.load_grouping_field(block, prim_group, "velocity", j,
				    velocity_j);
  array_factory.load_grouping_field(block, prim_group, "velocity", j,
				    velocity_k);
  array_factory.load_grouping_field(block, prim_group, "bfield", j, bfield_j);
  array_factory.load_grouping_field(block, prim_group, "bfield", k, bfield_k);

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

// Helper function to computes derives of x, y, and z with respect to r_{dim}.
// r_{dim} maps to x, y, and z for dim = 0, 1, 2
// The helper function computes dx/dr_{dim}, dy/dr_{dim}, dz/dr_{dim}
// We represent these with the shorthand dxddim, dyddim, dzddim
void aligned_dim_derivatives_(int dim, int &dxddim, int &dyddim, int &dzddim)
{
  if (dim == 0){
    dxddim = 1;
    dyddim = 0;
    dzddim = 0;
  } else if (dim == 1){
    dyddim = 1;
    dzddim = 0;
    dxddim = 0;
  } else {
    dzddim = 1;
    dxddim = 0;
    dyddim = 0;
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
  int dxdj, dydj, dzdj, dxdk, dydk, dzdk;
  aligned_dim_derivatives_((dim+1)%2, dxdj, dydj, dzdj);
  aligned_dim_derivatives_((dim+2)%2, dxdk, dydk, dzdk);

  // Initialize Cell-Centered E-fields
  EnzoArray<enzo_float> Ec, Ec_jp1, Ec_kp1, Ec_jkp1;
  EnzoFieldArrayFactory array_factory;
  array_factory.initialize_field_array(block, Ec, center_efield_id);
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
  array_factory.load_temp_interface_grouping_field(block, efield_group,
						   "efield", dim, Eedge,
						   dim == 0, dim == 1,
						   dim == 2);

  // Initialize face-centered E-fields
  EnzoArray<enzo_float> Ej, Ej_kp1, Ek, Ek_jp1;
  array_factory.load_temp_interface_grouping_field(block, jflux_group, "bfield",
						   dim, Ej, dxdj==0, dydj==0,
						   dzdj==0);
  Ej_kp1.initialize_subarray(Ej, dzdk, Ej.length_dim2(),
			         dydk, Ej.length_dim1(),
			         dxdk, Ej.length_dim0());
  array_factory.load_temp_interface_grouping_field(block, kflux_group, "bfield",
						   dim, Ek, dxdk==0, dydk==0,
						   dzdk==0);
  Ek_jp1.initialize_subarray(Ek, dzdj, Ek.length_dim2(),
			         dydj, Ek.length_dim1(),
			         dxdj, Ek.length_dim0());

  // Initialize the weight arrays
  EnzoArray<enzo_float> Wj, Wj_kp1, Wk, Wk_jp1;
  array_factory.load_temp_interface_grouping_field(block, weight_group,
						   "weight", (dim+1)%3, Wj,
						   dxdj==0, dydj==0, dzdj==0);
  Wj_kp1.initialize_subarray(Wj, dzdk, Wj.length_dim2(),
			         dydk, Wj.length_dim1(),
			         dxdk, Wj.length_dim0());
  array_factory.load_temp_interface_grouping_field(block, weight_group,
						   "weight", (dim+2)%3, Wk,
						   dxdk==0, dydk==0, dzdk==0);
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


// Compute the face-centered B-field component along the ith dimension
//
// Bnew_i(k, j, i+1/2) = Bold_i(k, j, i+1/2) -
//     dt/dj*(E_k(    k,j+1/2,i+1/2) - E_k(    k,j-1/2,i+1/2) +
//     dt/dk*(E_j(k+1/2,    j,i+1/2) - E_j(k-1/2,    j,i+1/2)
// [The positioning of dt/dj with respect to E_k is correct]
//
// Rewrite equation (to allow for handling of 2d and 3d grids):
// Bnew_i(k, j, i-1/2) =
//   Bold_i(k, j, i-1/2) - E_k_term(k,j,i+1/2) + E_j_term(k,j,i+1/2)
//
//   if (3D mesh || (2D mesh && dim==0)):
//     E_k_term(k,j,i+1/2) = dt/dj*(E_k(k,j+1/2,i+1/2) - E_k(k,j-1/2,i+1/2))
//   else:
//     E_k_term(k,j,i+1/2)  = 0
//
//   if (3D mesh || (2D mesh && dim==1)):
//     E_j_term(k,j,i+1/2) = dt/dk*(E_j(k+1/2,j,i+1/2) - E_j(k-1/2,j,i+1/2))
//   else:
//     E_j_term= 0
//
void EnzoConstrainedTransport::update_bfield(Block *block, int dim,
					     Grouping &efield_group,
					     Grouping &cur_bfieldi_group,
					     Grouping &out_bfieldi_group,
					     enzo_float dt)
{
  // determine alignment of j,k axes with respect to x,y, and z
  int  dxdj, dydj, dzdj, dxdk, dydk, dzdk;
  aligned_dim_derivatives_((dim+1)%3, dxdj, dydj, dzdj);
  aligned_dim_derivatives_((dim+2)%3, dxdk, dydk, dzdk);

  // determine if grid is 2D or 3D AND compute the ratios of dt to the 
  // widths of cells along j and k directions
  EnzoBlock *enzo_block = enzo::block(block);
  bool three_dim = enzo_block->GridDimension[2]>1;
  enzo_float dtdj = dt/enzo_block->CellWidth[(dim+1)%3];
  enzo_float dtdk = dt/enzo_block->CellWidth[(dim+2)%3];

  // Load interface bfields
  EnzoFieldArrayFactory array_factory;
  EnzoArray<enzo_float> cur_bfield, out_bfield;
  array_factory.load_interior_bfieldi_field(block, cur_bfieldi_group, dim,
					    cur_bfield);
  array_factory.load_interior_bfieldi_field(block, out_bfieldi_group, dim,
					    out_bfield);

  // Load edge centered efields
  EnzoArray<enzo_float> E_j, E_k;
  const bool use_E_j = (three_dim || dim == 1);
  // For 2d grid, only need to load E_k when dim == 1 (in that case E_j = E_z)
  if (use_E_j){
    array_factory.load_temp_interface_grouping_field(block, efield_group,
						     "efield", (dim+1)%3, E_j,
						     dxdj == 0, dydj == 0,
						     dzdj == 0);
  }
  const bool use_E_k = (three_dim || dim == 0);
  // For 2d grid, only need to load E_k when dim == 0 (in that case E_k = E_z)
  if (use_E_k){
    array_factory.load_temp_interface_grouping_field(block, efield_group,
						     "efield", (dim+2)%3, E_k,
						     dxdk == 0, dydk == 0,
						     dzdk == 0);
  }

  // Integration limits are compatible with 2D and 3D grids
  for (int iz=0; iz<out_bfield.length_dim2(); iz++) {
    for (int iy=0; iy<out_bfield.length_dim1(); iy++) {
      for (int ix=0; ix<out_bfield.length_dim0(); ix++) {

	enzo_float E_k_term, E_j_term;
	if (use_E_k){
	  // E_k_term(k,j,i+1/2) =
	  //   dt/dj*(E_k(k,j+1/2,i+1/2) - E_k(k,j-1/2,i+1/2))
	  E_k_term = dtdj*(E_k(iz,iy,ix) - E_k(iz-dzdj,iy-dydj,ix-dxdj));
	} else {
	  E_k_term = 0.;
	}

	if (use_E_j){
	  // E_j_term(k,j,i+1/2) =
	  //   dt/dk*(E_j(k+1/2,j,i+1/2) - E_j(k-1/2,j,i+1/2))
	  E_j_term = dtdk*(E_j(iz,iy,ix) - E_j(iz-dzdk,iy-dydk,ix-dxdk)); 
	} else {
	  E_j_term = 0.;
	}

	// Bnew_i(k, j, i-1/2) =
	//   Bold_i(k, j, i-1/2) - E_k_term(k,j,i+1/2) + E_j_term(k,j,i+1/2)
	out_bfield(iz,iy,ix) = cur_bfield(iz,iy,ix) - E_k_term + E_j_term;
      }
    }
  }
}


// Compute cell-centered bfield along dimension i
//   B_i(k,j,i) = 0.5*(B_i(k,j,i+1/2) + B_i(k,j,i-1/2))
// For a simpler implementation, we will rewrite this as:
//   B_i(k,j,i+1) = 0.5*(B_i(k,j,i+3/2) + B_i(k,j,i+1/2))
// We define:
//   B_center(k,j,i)   ->  B_i(k,j,i+1)
//   Bi_left(k,j,i)    ->  B_i(k,j,i+1/2)
//   Bi_right(k,j,i)   ->  B_i(k,j,i+3/2)
void EnzoConstrainedTransport::compute_center_bfield(Block *block, int dim,
						     Grouping &cons_group,
						     Grouping &bfieldi_group,
						     bool compute_outer)
{
  // Identify which axis dimension i is aligned with
  int  dxdi, dydi, dzdi;
  aligned_dim_derivatives_(dim, dxdi, dydi, dzdi);

  EnzoFieldArrayFactory array_factory;
  EnzoArray<enzo_float> b_center, bi_left, bi_right;
  int zstop, ystop, xstop;
  
  if (compute_outer){
    // In this case, compute the centered B-field for all cells in the grid
    // including cells at i = 0 and i = imax-1, where i is the index for the
    // axis aligned with dim.

    // Load cell-centerd field
    array_factory.load_grouping_field(block, cons_group, "bfield", dim,
				      b_center);

    // If we include the exterior interface b-fields, then we want to iterate
    // over the entire grid
    zstop = b_center.length_dim2();
    ystop = b_center.length_dim1();
    xstop = b_center.length_dim0();

    // Load Face-centered field - including values on the exterior faces
    array_factory.load_grouping_field(block, bfieldi_group, "bfield", dim,
				      bi_left);

  } else {
    // In this case, compute the centered B-field for all cells in the grid
    // except cells at i = 0 and i = imax-1, where i is the index for the
    // axis aligned with dim.
    
    // Load cell-centered B-field
    // b_center should not include any cells with ix=0, iy=0, or iz =0
    // since we do not have interface values at i-1/2 and don't need values
    // along the outmost layer
    EnzoArray<enzo_float> bfield;
    array_factory.load_grouping_field(block, cons_group, "bfield", dim, bfield);
    b_center.initialize_subarray(bfield, dzdi, bfield.length_dim2(),
				 dydi, bfield.length_dim1(),
				 dxdi, bfield.length_dim0());
    // Iteration Limits - need to stop when on the last index along dim i
    zstop = b_center.length_dim2() - dzdi;
    ystop = b_center.length_dim1() - dydi;
    xstop = b_center.length_dim0() - dxdi;
    
    // Load Face-centered field - excluding values on the exterior faces
    array_factory.load_interior_bfieldi_field(block, bfieldi_group, dim,
					      bi_left);
  }

  // Get the view of the Face-center field that starting from i=1
  bi_right.initialize_subarray(bi_left, dzdi, bi_left.length_dim2(),
			       dydi, bi_left.length_dim1(),
			       dxdi, bi_left.length_dim0());

  // iteration limits are compatible with a 2D grid and 3D grid
  for (int iz=0; iz<zstop; iz++) {
    for (int iy=0; iy<ystop; iy++) {
      for (int ix=0; ix<xstop; ix++) {

	b_center(iz,iy,ix) = 0.5*(bi_left(iz,iy,ix) + bi_right(iz,iy,ix));
      }
    }
  }
}
