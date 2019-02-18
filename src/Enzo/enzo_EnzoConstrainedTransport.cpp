#include "cello.hpp"
#include "enzo.hpp"

void EnzoConstrainedTransport::compute_center_efield (Block *block, int dim,
						      std::string center_efield_name,
						      Grouping &prim_group)
{
  // Load the E-field
  EnzoFieldArrayFactory array_factory(block);
  EFlt3DArray efield = array_factory.from_name(center_efield_name);

  EnzoPermutedCoordinates coord(dim);
  int j = coord.j_axis();
  int k = coord.k_axis();

  // Load the jth and kth components of the velocity and cell-centered bfield
  EFlt3DArray velocity_j, velocity_k, bfield_j, bfield_k;
  velocity_j = array_factory.from_grouping(prim_group, "velocity", j);
  velocity_k = array_factory.from_grouping(prim_group, "velocity", k);
  bfield_j = array_factory.from_grouping(prim_group, "bfield", j);
  bfield_k = array_factory.from_grouping(prim_group, "bfield", k);

  for (int iz=0; iz<efield.dim_size(2); iz++) {
    for (int iy=0; iy<efield.dim_size(1); iy++) {
      for (int ix=0; ix<efield.dim_size(0); ix++) {
	efield(iz,iy,ix) = (-velocity_j(iz,iy,ix) * bfield_k(iz,iy,ix) +
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
//  dEdj_r = dEdx(ix-1/4,iy-1/2) =
//       W_y(    ix,iy-1/2)  * (E(    ix,  iy-1) - E(ix-1/2,  iy-1)) +
//    (1-W_y(    ix,iy-1/2)) * (E(    ix,    iy) - E(ix-1/2,    iy))
//  dEdj_l = dEdx(ix-3/4,iy-1/2) =
//       W_y(  ix-1,iy-1/2)  * (E(ix-1/2,  iy-1) - E(  ix-1,  iy-1)) +
//    (1-W_y(  ix-1,iy-1/2)) * (E(ix-1/2,    iy) - E(  ix-1,    iy))
//  dEdk_r = dEdy(ix-1/2,iy-1/4) =
//       W_x(ix-1/2,   iy)  * (E(   ix-1,    iy) - E(  ix-1,iy-1/2)) +
//    (1-W_x(ix-1/2,   iy)) * (E(     ix,    iy) - E(    ix,iy-1/2))
//  dEdk_l = dEdy(ix-1/2,iy-3/4) =
//       W_x(ix-1/2, iy-1)  * (E(   ix-1,iy-1/2) - E(  ix-1,  iy-1)) +
//    (1-W_x(ix-1/2, iy-1)) * (E(     ix,iy-1/2) - E(    ix,  iy-1))
//
// This solution to the problem is easier to program if we reframe it as the
// calculation of Ez(ix+1/2, iy+1/2, iz). The above equations are then
// rewritten as:
//  dEdj_r = dEdx(ix+3/4,iy+1/2) =
//       W_y(  ix+1,iy+1/2)  * (E(  ix+1,    iy) - E(ix+1/2,    iy)) +
//    (1-W_y(  ix+1,iy+1/2)) * (E(  ix+1,  iy+1) - E(ix+1/2,  iy+1))
//  dEdj_l = dEdx(ix+1/4,iy+1/2) =
//       W_y(    ix,iy+1/2)  * (E(ix+1/2,    iy) - E(    ix,    iy)) +
//    (1-W_y(    ix,iy+1/2)) * (E(ix+1/2,  iy+1) - E(    ix,  iy+1))
//  dEdk_r = dEdy(ix+1/2,iy+3/4) =
//       W_x(ix+1/2,  iy+1)  * (E(     ix,  iy+1) - E(    ix,iy+1/2)) +
//    (1-W_x(ix+1/2,  iy+1)) * (E(   ix+1,  iy+1) - E(  ix+1,iy+1/2))
//  dEdk_l = dEdy(ix+1/2,iy+1/4) =
//       W_x(ix+1/2,   iy)  * (E(     ix,iy+1/2) - E(    ix,    iy)) +
//    (1-W_x(ix+1/2,   iy)) * (E(   ix+1,iy+1/2) - E(  ix+1,    iy))
//
// Now let's rewrite in totally Generalized notation:  [z->i, x->j, y->k] 
//  dEdj_r = dEdj(k+1/2,j+3/4,i) =
//       W_k(k+1/2,  j+1)  * (E(    k,  j+1) - E(    k,j+1/2)) +
//    (1-W_k(k+1/2,  j+1)) * (E(  k+1,  j+1) - E(  k+1,j+1/2))
//  dEdj_l = dEdj(k+1/2,j+1/4,i) =
//       W_k(k+1/2,    j)  * (E(    k,j+1/2) - E(    k,    j)) +
//    (1-W_k(k+1/2,    j)) * (E(  k+1,j+1/2) - E(  k+1,    j))
//  dEdk_r = dEdk(k+3/4,j+1/2,i) =
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
		   EFlt3DArray &Eedge, EFlt3DArray &Wj, EFlt3DArray &Wj_kp1,
		   EFlt3DArray &Wk, EFlt3DArray &Wk_jp1, EFlt3DArray &Ec,
		   EFlt3DArray &Ec_jkp1, EFlt3DArray &Ec_jp1,
		   EFlt3DArray &Ec_kp1, EFlt3DArray &Ej, EFlt3DArray &Ej_kp1,
		   EFlt3DArray &Ek, EFlt3DArray &Ek_jp1)
{
  for (int iz = zstart; iz < zstop; iz++){
    for (int iy = ystart; iy < ystop; iy++){
      for (int ix = xstart; ix < xstop; ix++){

	enzo_float dEdj_l, dEdj_r, dEdk_l, dEdk_r;

	//  dEdj(k+1/2,j+3/4,i) =
	//       W_k(k+1/2,  j+1)  * (E(    k,  j+1) - E(    k,j+1/2)) +
	//    (1-W_k(k+1/2,  j+1)) * (E(  k+1,  j+1) - E(  k+1,j+1/2))
	dEdj_r =
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
	dEdk_r =
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
				(dEdj_l-dEdj_r) + (dEdk_l - dEdk_r));

	//if (iy == 4 && iz == 4 && ix == 4){
	//  CkPrintf("Wj = %.15g, Wk = %.15g\n",
	//	    Wj(iz,iy,ix), Wk(iz,iy,ix));
	//  CkPrintf(("dEdk_l = %.15g, dEdk_r = %.15g\n"
	//	    "dEdj_l = %.15g, dEdj_r = %.15g\n"),
	//	    dEdk_l, dEdk_r, dEdj_l,dEdj_r);
	//  CkPrintf(("Ek = %.15g, Ek_jp1 = %.15g\n"
	//	    "Ej = %.15g, Ej_kp1 = %.15g\n"),
	//	   Ek(iz,iy,ix), Ek_jp1(iz,iy,ix),
	//	   Ej(iz,iy,ix),Ej_kp1(iz,iy,ix));
	//  CkPrintf("Eedge = %.15g\n", Eedge(iz,iy,ix));
	//  fflush(stdout);
	//}
      }
    }
  }
}


void inplace_entry_multiply_(EFlt3DArray &array, enzo_float val){
  for (int iz = 0; iz < array.dim_size(2); iz++){
    for (int iy = 0; iy < array.dim_size(1); iy++){
      for (int ix = 0; ix < array.dim_size(0); ix++){
	array(iz,iy,ix)*=val;
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
						    std::string center_efield_name,
						    Grouping &efield_group,
						    Grouping &jflux_group,
						    Grouping &kflux_group,
						    Grouping &prim_group,
						    Grouping &weight_group)
{

  EnzoPermutedCoordinates coord(dim);
  // determine components of j and k unit vectors:
  int j_x, j_y, j_z, k_x, k_y, k_z;
  coord.j_unit_vector(j_x, j_y, j_z);
  coord.k_unit_vector(k_x, k_y, k_z);

  // Initialize Cell-Centered E-fields
  EFlt3DArray Ec, Ec_jp1, Ec_kp1, Ec_jkp1;
  EnzoFieldArrayFactory array_factory(block);
  Ec = array_factory.from_name(center_efield_name);
  Ec_jp1  = coord.left_edge_offset(Ec, 0, 1, 0);
  Ec_kp1  = coord.left_edge_offset(Ec, 1, 0, 0);
  Ec_jkp1 = coord.left_edge_offset(Ec, 1, 1, 0);

  // Initialize edge-centered Efield [it maps (k,j,i) -> (k+1/2,j+1/2,i)]
  EFlt3DArray Eedge = array_factory.from_grouping(efield_group, "efield", dim);

  // Initialize face-centered E-fields
  EFlt3DArray Ej, Ej_kp1, Ek, Ek_jp1;

  // Ex(k,j+1/2,i) = -1.*yflux(Bz)
  Ej = array_factory.from_grouping(jflux_group, "bfield", coord.k_axis());
  inplace_entry_multiply_(Ej,-1.);
  Ej_kp1 = coord.left_edge_offset(Ej, 1, 0, 0);

  // Ex(k+1/2,j,i) = zflux(By)
  Ek = array_factory.from_grouping(kflux_group, "bfield", coord.j_axis());
  // No need to multiply entries by minus 1
  Ek_jp1 = coord.left_edge_offset(Ek, 0, 1, 0);

  // Initialize the weight arrays
  EFlt3DArray Wj, Wj_kp1, Wk, Wk_jp1;
  Wj = array_factory.from_grouping(weight_group, "weight", coord.j_axis());
  Wj_kp1 = coord.left_edge_offset(Wj, 1, 0, 0);
  Wk = array_factory.from_grouping(weight_group, "weight", coord.k_axis());
  Wk_jp1 = coord.left_edge_offset(Wk, 0, 1, 0);

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
  //    2D grid with i aligned along z:  istart = 0     istop = imax - 1
  //    3D grid:                         istart = 1     istop = imax - 1
  //    In all cases:                    jstart = 0     jstop = jmax - 1
  //                                     kstart = 0     kstop = kmax - 1

  int xstart = 1 - j_x - k_x; // if dim==0: 1, otherwise: 0
  int ystart = 1 - j_y - k_y;
  int zstart = (Ec.dim_size(2) == 1) ? 0 : 1 - j_z - k_z;

  int zstop = (Ec.dim_size(2) == 1) ? 1 : Ec.dim_size(2) - 1;

  compute_edge_(xstart, ystart, zstart,
		Ec.dim_size(0) - 1, Ec.dim_size(1) - 1, zstop,
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
  EnzoPermutedCoordinates coord(dim);
  // determine components of j and k unit vectors:
  int j_x, j_y, j_z, k_x, k_y, k_z;
  coord.j_unit_vector(j_x, j_y, j_z);
  coord.k_unit_vector(k_x, k_y, k_z);

  // determine if grid is 2D or 3D AND compute the ratios of dt to the 
  // widths of cells along j and k directions
  EnzoBlock *enzo_block = enzo::block(block);
  bool three_dim = enzo_block->GridDimension[2]>1;
  enzo_float dtdj = dt/enzo_block->CellWidth[coord.j_axis()];
  enzo_float dtdk = dt/enzo_block->CellWidth[coord.k_axis()];

  // Load interface bfields
  EnzoFieldArrayFactory array_factory(block);
  EFlt3DArray cur_bfield, out_bfield;
  cur_bfield = array_factory.interior_bfieldi(cur_bfieldi_group, dim);
  out_bfield = array_factory.interior_bfieldi(out_bfieldi_group, dim);

  // Load edge centered efields
  EFlt3DArray E_j, E_k;
  // For 2D grid: only need to load E_k when dim == 1 (in that case E_j = E_z)
  //              only need to load E_k when dim == 0 (in that case E_k = E_z)
  const bool use_E_j = (three_dim || dim == 1);
  const bool use_E_k = (three_dim || dim == 0);
  if (use_E_j){
    E_j = array_factory.from_grouping(efield_group, "efield", coord.j_axis());
  }
  if (use_E_k){
    E_k = array_factory.from_grouping(efield_group, "efield", coord.k_axis());
  }

  int zstart = 0;
  int zstop = 1;
  if (three_dim){
    zstart = j_z + k_z;
    zstop = out_bfield.dim_size(2)-j_z-k_z;
  }

  // We could simplify this iteration by using subarrays - However, it would be
  // more complicated
  for (int iz=zstart; iz<zstop; iz++) {
    for (int iy=j_y+k_y; iy<out_bfield.dim_size(1)-j_y-k_y; iy++) {
      for (int ix=j_x+k_x; ix<out_bfield.dim_size(0)-j_x-k_x; ix++) {

	enzo_float E_k_term, E_j_term;
	if (use_E_k){
	  // E_k_term(k,j,i+1/2) =
	  //   dt/dj*(E_k(k,j+1/2,i+1/2) - E_k(k,j-1/2,i+1/2))
	  E_k_term = dtdj*(E_k(iz,iy,ix) - E_k(iz-j_z,iy-j_y,ix-j_x));
	} else {
	  E_k_term = 0.;
	}

	if (use_E_j){
	  // E_j_term(k,j,i+1/2) =
	  //   dt/dk*(E_j(k+1/2,j,i+1/2) - E_j(k-1/2,j,i+1/2))
	  E_j_term = dtdk*(E_j(iz,iy,ix) - E_j(iz-k_z,iy-k_y,ix-k_x)); 
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

// This method also intentionally includes calculation of bfields in the
// outermost cells so that it can be used to initially setup the bfield.
//
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
						     Grouping &bfieldi_group)
{
  EnzoPermutedCoordinates coord(dim);
  EnzoFieldArrayFactory array_factory(block);

  // Load cell-centerd field
  EFlt3DArray b_center = array_factory.from_grouping(cons_group, "bfield",
						     coord.i_axis());
  // Load Face-centered fields
  EFlt3DArray bi_left = array_factory.from_grouping(bfieldi_group, "bfield",
						    coord.i_axis());
  // Get the view of the Face-center field that starting from i=1
  EFlt3DArray bi_right = coord.left_edge_offset(bi_left,0,0,1);

  // iteration limits are compatible with a 2D grid and 3D grid
  for (int iz=0; iz<b_center.dim_size(2); iz++) {
    for (int iy=0; iy<b_center.dim_size(1); iy++) {
      for (int ix=0; ix<b_center.dim_size(0); ix++) {
	b_center(iz,iy,ix) = 0.5*(bi_left(iz,iy,ix) + bi_right(iz,iy,ix));
      }
    }
  }
}
