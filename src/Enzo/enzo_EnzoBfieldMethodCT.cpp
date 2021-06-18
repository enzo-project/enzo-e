// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBfieldMethodCT.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon May 6 2019
/// @brief    [\ref Enzo] Implementation of EnzoBfieldMethodCT

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoBfieldMethodCT::register_target_block_
(Block *block, bool first_initialization) noexcept
{
  // setup bfieldi_l_ (initialize arrays that wrap the fields holding each
  // component of the interface centered magnetic field).
  const std::string field_names[] = {"bfieldi_x", "bfieldi_y", "bfieldi_z"};
  EnzoFieldArrayFactory array_factory(block, 0); // stale_depth = 0
  for (std::size_t i = 0; i<3; i++){
    bfieldi_l_[i] = array_factory.from_name(field_names[i]);
  }

  EnzoBlock *enzo_block = enzo::block(block);
  cell_widths_ = (const enzo_float*)enzo_block->CellWidth;

  if (first_initialization){
    // if num_partial_timesteps() > 1, setup temp_bfieldi_l_ (they have the same
    // shape as each component of bfieldi_l_ and manage their own memory).
    if (num_partial_timesteps() > 1){
      for (std::size_t i = 0; i<3; i++){
        temp_bfieldi_l_[i] = EFlt3DArray(bfieldi_l_[i].shape(0),
                                         bfieldi_l_[i].shape(1),
                                         bfieldi_l_[i].shape(2));
      }
    }

    // get the standard size of cell-centered fields
    int mz = bfieldi_l_[0].shape(0); // interface bfield_z cell-centered along z
    int my = bfieldi_l_[0].shape(1); // interface bfield_y cell-centered along y
    // interface bfield_x is face-centered along x (and the array includes
    // values on external faces of the grid)
    int mx = bfieldi_l_[0].shape(2) - 1;

    // allocate arrays to hold weight values for each dimension. For a given
    // dimension, the array tracks the upwind/downwind direction on the cell
    // interfaces for that dimension (they exclude exterior faces of the block)
    weight_l_[0] = EFlt3DArray(  mz,   my, mx-1);
    weight_l_[1] = EFlt3DArray(  mz, my-1,   mx);
    weight_l_[2] = EFlt3DArray(mz-1,   my,   mx);

    // allocate arrays to hold for each component of the edge-centered electric
    // fields. The array for component i is cell-centered along i and face
    // centered along j and k.
    edge_efield_l_[0] = EFlt3DArray(mz-1, my-1,   mx);
    edge_efield_l_[1] = EFlt3DArray(mz-1,   my, mx-1);
    edge_efield_l_[2] = EFlt3DArray(  mz, my-1, mx-1);

    // allocate array to temporarily store cell-centered e-field
    center_efield_ = EFlt3DArray(mz,my,mx);
  }
}

//----------------------------------------------------------------------

void EnzoBfieldMethodCT::check_required_fields() const noexcept
{
  FieldDescr * field_descr = cello::field_descr();
  ASSERT("EnzoBfieldMethodCT::check_required_fields",
	 ("There must be face-centered permanent fields called \"bfieldi_x\","
	  "\"bfield_y\" and \"bfield_z\"."),
	 (field_descr->is_field("bfield_x") &&
	  field_descr->is_field("bfield_y") &&
	  field_descr->is_field("bfield_z")));

  std::vector<std::string> names = {"bfieldi_x", "bfieldi_y", "bfieldi_z"};
  std::vector<std::string> axes = {"x", "y", "z"};
  for (std::size_t i = 0; i < 3; i++){
    std::string name = names[i];
    // first check that field exists
    ASSERT1("EnzoBfieldMethodCT::check_required_fields",
	    "There must be face-centered permanent fields called \"%s\"",
	    name.c_str(), field_descr->is_field(name));

    int field_id = field_descr->field_id(name);
    // next check the centering of the field
    int centering[3] = {0, 0, 0};
    field_descr->centering(field_id, &centering[0], &centering[1],
			   &centering[2]);
    for (std::size_t j = 0; j<3; j++){
      if (j!=i){
	ASSERT2("EnzoBfieldMethodCT::update_refresh",
		"The \"%s\" field must be cell-centered along the %s-axis",
		name.c_str(), axes[j].c_str(), centering[j] == 0);
      } else {
	ASSERT2("EnzoBfieldMethodCT::update_refresh",
		"The \"%s\" field must be face-centered along the %s-axis",
		name.c_str(), axes[j].c_str(), centering[j] == 1);
      }
    }

  }
}

//----------------------------------------------------------------------

void EnzoBfieldMethodCT::correct_reconstructed_bfield
(EnzoEFltArrayMap &l_map, EnzoEFltArrayMap &r_map, int dim,
 int stale_depth) noexcept
{
  require_registered_block_(); // confirm that target_block_ is valid

  std::array<EFlt3DArray,3> *cur_bfieldi_l;
  if (partial_timestep_index() == 0){
    cur_bfieldi_l = &bfieldi_l_;
  } else if (partial_timestep_index() == 1){
    cur_bfieldi_l = &temp_bfieldi_l_;
  } else {
    ERROR1("EnzoBfieldMethodCT::correct_reconstructed_bfield",
	   "Unsure how to handle current partial_timestep_index_ of %d",
	   partial_timestep_index());
  }

  if ((dim < 0) || (dim > 2)){
    ERROR("EnzoBfieldMethodCT::correct_reconstructed_bfield",
          "dim has an invalid value");
  } else {
    // the interface bfield values held in *cur_bfieldi_l[dim] includes values
    // on the exterior faces of the block. We only need the values on the
    // interior faces.
    EnzoPermutedCoordinates coord(dim);
    CSlice full_ax(nullptr,nullptr);
    EFlt3DArray bfield = coord.get_subarray((*cur_bfieldi_l)[dim],
                                            full_ax, full_ax, CSlice(1,-1));

    const std::string names[3] = {"bfield_x", "bfield_y", "bfield_z"};
    EFlt3DArray l_bfield = l_map.at(names[dim]);
    EFlt3DArray r_bfield = r_map.at(names[dim]);

    // All 3 array objects are the same shape
    for (int iz = stale_depth; iz< bfield.shape(0) - stale_depth; iz++) {
      for (int iy = stale_depth; iy< bfield.shape(1) - stale_depth; iy++) {
        for (int ix = stale_depth; ix < bfield.shape(2) - stale_depth; ix++) {
          l_bfield(iz,iy,ix) = bfield(iz,iy,ix);
          r_bfield(iz,iy,ix) = bfield(iz,iy,ix);
        }
      }
    }

  }
}

//----------------------------------------------------------------------

void EnzoBfieldMethodCT::identify_upwind(const EnzoEFltArrayMap &flux_map,
                                         int dim, int stale_depth) noexcept
{
  require_registered_block_(); // confirm that target_block_ is valid

  //  - Currently, weight is set to 1.0 if upwind is in positive direction of
  //    the current dimension, 0 if upwind is in the negative direction of the
  //    current dimension, or 0.5 if there is no upwind direction
  //  - At present, the weights are unnecessary (the same information is
  //    encoded in density flux to figure out this information). However, this
  //    functionallity is implemented in case we decide to adopt the weighting
  //    scheme from Athena++, which requires knowledge of the reconstructed
  //    densities.

  if ((dim < 0) || (dim > 2)){
    ERROR("EnzoBfieldMethodCT::identify_upwind",
          "dim has an invalid value");
  } else {
    EFlt3DArray density_flux = flux_map.get("density", stale_depth);

    CSlice stale_slc = (stale_depth > 0) ?
      CSlice(stale_depth,-stale_depth) : CSlice(nullptr, nullptr);

    EFlt3DArray weight_field = weight_l_[dim].subarray(stale_slc, stale_slc,
                                                       stale_slc);

    // Iteration limits compatible with both 2D and 3D grids
    for (int iz=0; iz<density_flux.shape(0); iz++) {
      for (int iy=0; iy<density_flux.shape(1); iy++) {
        for (int ix=0; ix<density_flux.shape(2); ix++) {
          // density flux is the face-centered density times the face-centered
          // velocity along dim

          if ( density_flux(iz,iy,ix) > 0){
            weight_field(iz,iy,ix) = 1.0;
          } else if ( density_flux(iz,iy,ix) < 0){
            weight_field(iz,iy,ix) = 0.0;
          } else {
            weight_field(iz,iy,ix) = 0.5;
          }
        }
      }
    }

  }
}

//----------------------------------------------------------------------

void EnzoBfieldMethodCT::update_all_bfield_components
(EnzoEFltArrayMap &cur_prim_map, EnzoEFltArrayMap &xflux_map,
 EnzoEFltArrayMap &yflux_map, EnzoEFltArrayMap &zflux_map,
 EnzoEFltArrayMap &out_centered_bfield_map, enzo_float dt,
 int stale_depth) noexcept
{
  // The current interface values of the interface magnetic fields that are to
  // be updated are always the values from the start of the timestep
  std::array<EFlt3DArray,3> *cur_bfieldi_l = &bfieldi_l_;
  // Now determine the set of arrays where the updated interface bfields are to
  // be stored:
  std::array<EFlt3DArray,3> *out_bfieldi_l;
  if ((partial_timestep_index()== 1) || (num_partial_timesteps() == 1)){
    out_bfieldi_l = &bfieldi_l_; // the same as cur_bfieldi_l
  } else if (partial_timestep_index()== 0){
    out_bfieldi_l = &temp_bfieldi_l_;
  } else {
    ERROR1("EnzoBfieldMethodCT::correct_reconstructed_bfield",
	   "Unsure how to handle current partial_timestep_index_ of %d",
	   partial_timestep_index());
  }

  // First, compute the edge-centered Electric fields (each time, it uses
  // the current integrable quantities)
  EnzoBfieldMethodCT::compute_all_edge_efields(cur_prim_map, xflux_map,
                                                     yflux_map, zflux_map,
                                                     center_efield_,
                                                     edge_efield_l_, weight_l_,
                                                     stale_depth);

  // Update longitudinal B-field (add source terms of constrained transport)
  for (int dim = 0; dim<3; dim++){
    EnzoBfieldMethodCT::update_bfield
      (cell_widths_, dim, edge_efield_l_, (*cur_bfieldi_l)[dim],
       (*out_bfieldi_l)[dim], dt, stale_depth);
  }

  // Finally, update cell-centered B-field
  const std::string names[3] = {"bfield_x", "bfield_y", "bfield_z"};
  for (int dim = 0; dim<3; dim++){
    EnzoBfieldMethodCT::compute_center_bfield
      (dim, out_centered_bfield_map[names[dim]], (*out_bfieldi_l)[dim],
       stale_depth);
  }
}

//----------------------------------------------------------------------

void EnzoBfieldMethodCT::compute_center_efield
(int dim, EFlt3DArray &efield, const EnzoEFltArrayMap &prim_map,
 int stale_depth)
{
  const std::string v_names[3] = {"velocity_x", "velocity_y", "velocity_z"};
  const std::string b_names[3] = {"bfield_x", "bfield_y", "bfield_z"};

  EnzoPermutedCoordinates coord(dim);
  int j = coord.j_axis();
  int k = coord.k_axis();

  // Load the jth and kth components of the velocity and cell-centered bfield
  EFlt3DArray velocity_j = prim_map.at(v_names[j]);
  EFlt3DArray velocity_k = prim_map.at(v_names[k]);
  EFlt3DArray bfield_j = prim_map.at(b_names[j]);
  EFlt3DArray bfield_k = prim_map.at(b_names[k]);

  for (int iz=stale_depth; iz<efield.shape(0)-stale_depth; iz++) {
    for (int iy=stale_depth; iy<efield.shape(1)-stale_depth; iy++) {
      for (int ix=stale_depth; ix<efield.shape(2)-stale_depth; ix++) {
	efield(iz,iy,ix) = (-velocity_j(iz,iy,ix) * bfield_k(iz,iy,ix) +
                            velocity_k(iz,iy,ix) * bfield_j(iz,iy,ix));
      }
    }
  }
}

//----------------------------------------------------------------------

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
//     Ec_jp1(k,j,i)   =     E(    k,  j+1,i)
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

      }
    }
  }
}

//----------------------------------------------------------------------

void inplace_entry_multiply_(EFlt3DArray &array, enzo_float val){
  for (int iz = 0; iz < array.shape(0); iz++){
    for (int iy = 0; iy < array.shape(1); iy++){
      for (int ix = 0; ix < array.shape(2); ix++){
	array(iz,iy,ix)*=val;
      }
    }
  }
}

//----------------------------------------------------------------------

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
// The Athena++ code calculates a quantity they refer to as v_over_c at all
// cell-faces.
//   - Basically, this tells them how much to weight derivatives while
//     computing the edge_efield. If deemed necessary, these values can be
//     passed as weight_l
//   - Weight_l includes 3 temporary arrays centered on the faces looking
//     down the x, y, and z direction. Currently, it expects values of 1 and 0
//     to indicate that the upwind direction is in positive and negative
//     direction, or 0 to indicate no upwind direction.
void EnzoBfieldMethodCT::compute_edge_efield
(int dim, EFlt3DArray &center_efield, EFlt3DArray &edge_efield,
 EnzoEFltArrayMap &jflux_map, EnzoEFltArrayMap &kflux_map,
 std::array<EFlt3DArray,3> &weight_l, int stale_depth)
{

  EnzoPermutedCoordinates coord(dim);
  // determine components of j and k unit vectors:
  int j_x, j_y, j_z, k_x, k_y, k_z;
  coord.j_unit_vector(j_z, j_y, j_x);
  coord.k_unit_vector(k_z, k_y, k_x);

  CSlice stale_slc = (stale_depth > 0) ?
      CSlice(stale_depth,-stale_depth) : CSlice(nullptr, nullptr);

  // Initialize Cell-Centered E-fields
  EFlt3DArray Ec, Ec_jp1, Ec_kp1, Ec_jkp1;
  Ec = center_efield.subarray(stale_slc, stale_slc, stale_slc);
  Ec_jp1  = coord.left_edge_offset(Ec, 0, 1, 0);
  Ec_kp1  = coord.left_edge_offset(Ec, 1, 0, 0);
  Ec_jkp1 = coord.left_edge_offset(Ec, 1, 1, 0);

  // Initialize edge-centered Efield [it maps (k,j,i) -> (k+1/2,j+1/2,i)]
  EFlt3DArray Eedge = edge_efield.subarray(stale_slc, stale_slc, stale_slc);

  // Initialize face-centered E-fields
  const std::string keys[] = {"bfield_x", "bfield_y", "bfield_z"};

  // Ex(k,j+1/2,i) = -1.*yflux(Bz)
  EFlt3DArray Ej = jflux_map.get(keys[coord.k_axis()],stale_depth);
  inplace_entry_multiply_(Ej,-1.);
  EFlt3DArray Ej_kp1 = coord.left_edge_offset(Ej, 1, 0, 0);

  // Ex(k+1/2,j,i) = zflux(By)
  EFlt3DArray Ek = kflux_map.get(keys[coord.j_axis()],stale_depth);
  // No need to multiply entries by minus 1
  EFlt3DArray Ek_jp1 = coord.left_edge_offset(Ek, 0, 1, 0);

  // Initialize the weight arrays
  EFlt3DArray Wj, Wj_kp1, Wk, Wk_jp1;
  Wj = weight_l[coord.j_axis()].subarray(stale_slc, stale_slc, stale_slc);
  Wj_kp1 = coord.left_edge_offset(Wj, 1, 0, 0);
  Wk = weight_l[coord.k_axis()].subarray(stale_slc, stale_slc, stale_slc);
  Wk_jp1 = coord.left_edge_offset(Wk, 0, 1, 0);

  // Integration limits
  //
  // If computing the edge E-field along z-direction:
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
  //    istart = 1     istop = imax - 1
  //    jstart = 0     jstop = jmax - 1
  //    kstart = 0     kstop = kmax - 1

  int xstart = 1 - j_x - k_x; // if dim==0: 1, otherwise: 0
  int ystart = 1 - j_y - k_y; // if dim==1: 1, otherwise: 0
  int zstart = 1 - j_z - k_z; // if dim==2: 1, otherwise: 0

  compute_edge_(xstart, ystart, zstart,
		Ec.shape(2) - 1, Ec.shape(1) - 1, Ec.shape(0) - 1,
		Eedge, Wj, Wj_kp1, Wk, Wk_jp1, Ec, Ec_jkp1, Ec_jp1, Ec_kp1,
		Ej, Ej_kp1, Ek, Ek_jp1);
}

//----------------------------------------------------------------------

void EnzoBfieldMethodCT::compute_all_edge_efields
  (EnzoEFltArrayMap &prim_map, EnzoEFltArrayMap &xflux_map,
   EnzoEFltArrayMap &yflux_map, EnzoEFltArrayMap &zflux_map,
   EFlt3DArray &center_efield, std::array<EFlt3DArray,3> &edge_efield_l,
   std::array<EFlt3DArray,3> &weight_l, int stale_depth)
{

  for (int i = 0; i < 3; i++){
    EnzoBfieldMethodCT::compute_center_efield(i, center_efield, prim_map,
                                                    stale_depth);

    EnzoEFltArrayMap *jflux_map;
    EnzoEFltArrayMap *kflux_map;
    if (i == 0){
      jflux_map = &yflux_map;
      kflux_map = &zflux_map;
    } else if (i==1){
      jflux_map = &zflux_map;
      kflux_map = &xflux_map;
    } else {
      jflux_map = &xflux_map;
      kflux_map = &yflux_map;
    }

    EnzoBfieldMethodCT::compute_edge_efield(i, center_efield,
                                                  edge_efield_l[i], *jflux_map,
                                                  *kflux_map, weight_l,
                                                  stale_depth);
  }
}

//----------------------------------------------------------------------

// Compute the face-centered B-field component along the ith dimension
//
// Bnew_i(k, j, i-1/2) = Bold_i(k, j, i-1/2) -
//     dt/dj*(E_k(    k,j+1/2,i-1/2) - E_k(    k,j-1/2,i-1/2) +
//     dt/dk*(E_j(k+1/2,    j,i-1/2) - E_j(k-1/2,    j,i-1/2)
// [The positioning of dt/dj with respect to E_k is correct]
//
// Bnew_i(k, j, i+1/2) =
//   Bold_i(k, j, i+1/2) - E_k_term(k,j,i+1/2) + E_j_term(k,j,i+1/2)
//
// Assuming 3D:
//   E_k_term(k,j,i+1/2) = dt/dj*(E_k(k,j+1/2,i+1/2) - E_k(k,j-1/2,i+1/2))
//   E_j_term(k,j,i+1/2) = dt/dk*(E_j(k+1/2,j,i+1/2) - E_j(k-1/2,j,i+1/2))
// 
// We define: (notation DIFFERENT from compute_edge_efield & compute_edge_)
//   ej_Lk(k,j,i) = E_j(k-1/2,     j, i+1/2) 
//   ej_Rk(k,j,i) = E_j(k+1/2,     j, i+1/2)
//   ek_Lj(k,j,i) = E_k(    k, j-1/2, i+1/2) 
//   ek_Rj(k,j,i) = E_k(    k, j+1/2, i+1/2)
//
// Then:
//   E_k_term(k,j,i+1/2) = dt/dj*(ek_Rj(k,j,i) - ek_Lj(k,j,i))
//   E_j_term(k,j,i+1/2) = dt/dk*(ej_Rk(k,j,i) - ej_Lk(k,j,i))
void EnzoBfieldMethodCT::update_bfield
(const enzo_float* &cell_widths, int dim,
 const std::array<EFlt3DArray,3> &efield_l,
 EFlt3DArray &cur_interface_bfield, EFlt3DArray &out_interface_bfield,
 enzo_float dt, int stale_depth)
{
  EnzoPermutedCoordinates coord(dim);

  // compute the ratios of dt to the widths of cells along j and k directions
  enzo_float dtdj = dt/cell_widths[coord.j_axis()];
  enzo_float dtdk = dt/cell_widths[coord.k_axis()];

  CSlice stale_slc = (stale_depth > 0) ?
      CSlice(stale_depth,-stale_depth) : CSlice(nullptr, nullptr);

  // Load edge centered efields 
  EFlt3DArray E_j, ej_Lk, ej_Rk, E_k, ek_Lj, ek_Rj;
  E_j = efield_l[coord.j_axis()].subarray(stale_slc, stale_slc, stale_slc);
  E_k = efield_l[coord.k_axis()].subarray(stale_slc, stale_slc, stale_slc);

  // Load interface bfields (includes exterior faces)
  EFlt3DArray cur_bfield, bcur, out_bfield, bout;
  cur_bfield = cur_interface_bfield.subarray(stale_slc, stale_slc, stale_slc);
  out_bfield = out_interface_bfield.subarray(stale_slc, stale_slc, stale_slc);

  // The following all assume that we are talking about the unstaled region
  // (i.e. we have dropped all staled cells)

  // Now to take slices. If the unstaled region of the grid has shape
  // (mk, mj, mi) then:
  //   - E_j has shape (mk-1, mj, mi-1)
  //       ej_Lk includes k=1/2 up to (but not including) k=mk-3/2
  //       ej_Rk includes k=3/2 up to (but not including) k=mk-1/2
  //   - E_k has shape (mk, mj-1, mi-1)
  //       ek_Lj includes j=1/2 up to (but not including) j=mj-3/2
  //       ek_Rj includes j=3/2 up to (but not including) j=mj-1/2
  //   - cur_bfield and out_bfield each have shape (mk, mj, mi+1)
  //       bnew and bout only include interior faces
  //
  // Also need to omit outermost layer of cell-centered vals
  CSlice full_ax(nullptr, nullptr); // includes full axis
  CSlice inner_cent(1,-1);          // excludes outermost cell-centered values

  // the following arrays should all have the same shape
  ej_Lk = coord.get_subarray(E_j, CSlice(0,      -1), inner_cent, full_ax);
  ej_Rk = coord.get_subarray(E_j, CSlice(1, nullptr), inner_cent, full_ax);
  ek_Lj = coord.get_subarray(E_k, inner_cent, CSlice(0,      -1), full_ax);
  ek_Rj = coord.get_subarray(E_k, inner_cent, CSlice(1, nullptr), full_ax);
  bcur = coord.get_subarray(cur_bfield, inner_cent, inner_cent, CSlice(1,-1));
  bout = coord.get_subarray(out_bfield, inner_cent, inner_cent, CSlice(1,-1));

  // We could simplify this iteration by using subarrays - However, it would be
  // more complicated
  for (int iz=0; iz<bout.shape(0); iz++) {
    for (int iy=0; iy<bout.shape(1); iy++) {
      for (int ix=0; ix<bout.shape(2); ix++) {

	// E_k_term(k,j,i+1/2) = dt/dj*(E_k(k,j+1/2,i+1/2) - E_k(k,j-1/2,i+1/2))
	enzo_float E_k_term = dtdj*(ek_Rj(iz,iy,ix) - ek_Lj(iz,iy,ix));

	// E_j_term(k,j,i+1/2) = dt/dk*(E_j(k+1/2,j,i+1/2) - E_j(k-1/2,j,i+1/2))
	enzo_float E_j_term = dtdk*(ej_Rk(iz,iy,ix) - ej_Lk(iz,iy,ix)); 

	// Bnew_i(k, j, i+1/2) =
	//   Bold_i(k, j, i+1/2) - E_k_term(k,j,i+1/2) + E_j_term(k,j,i+1/2)
	bout(iz,iy,ix) = bcur(iz,iy,ix) - E_k_term + E_j_term;
      }
    }
  }
}

//----------------------------------------------------------------------

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
void EnzoBfieldMethodCT::compute_center_bfield(int dim,
                                               EFlt3DArray &bfieldc_comp,
                                               EFlt3DArray &bfieldi_comp,
                                               int stale_depth)
{
  EnzoPermutedCoordinates coord(dim);
  CSlice stale_slc = (stale_depth > 0) ?
      CSlice(stale_depth,-stale_depth) : CSlice(nullptr, nullptr);

  // Load Cell-centered fields
  EFlt3DArray b_center = bfieldc_comp.subarray(stale_slc,stale_slc,stale_slc);
  // Load Face-centered fields
  EFlt3DArray bi_left = bfieldi_comp.subarray(stale_slc,stale_slc,stale_slc);
  // Get the view of the Face-center B-field that starting from i=1
  EFlt3DArray bi_right = coord.left_edge_offset(bi_left,0,0,1);

  // iteration limits are compatible with a 2D grid and 3D grid
  for (int iz=0; iz<b_center.shape(0); iz++) {
    for (int iy=0; iy<b_center.shape(1); iy++) {
      for (int ix=0; ix<b_center.shape(2); ix++) {
	b_center(iz,iy,ix) = 0.5*(bi_left(iz,iy,ix) + bi_right(iz,iy,ix));
      }
    }
  }
}
