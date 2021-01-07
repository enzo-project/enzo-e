// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoConstrainedTransport.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon May 6 2019
/// @brief    [\ref Enzo] Implementation of EnzoConstrainedTransport

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

// Helper function used to prepare groupings of temporary fields. These
// groupings only contain 3 fields - corresponding to vector components
//  - These vectors are all expected to be face-centered or edge centered.
//    (with the direction of centering depending on the dimension of the
//     component)
//  - face-centered: [face_center = True] fields are all centered on the face
//    corresponding to the direction of the vector component (e.g. x component
//    of weights is centered on the faces along the x dimension)
//  - edge-centered: [face_center = False] fields are centered on the edge
//    corresponding to the direction not pointed to by the vector component
//    (e.g. z-component of edge E-field is centered on the x and y edges)
//  - If exterior_faces is true, then the field will include space for the
//    exterior faces. If it is false, then the field will not include space for
//    the exterior faces.
//  - the names of the temporary fields are given by field_prefix + x,
//    field_prefix + y and field_prefix + z

void prep_temp_vector_grouping_(Field &field, std::string group_name,
				Grouping &grouping, std::string field_prefix,
				bool face_center, bool exterior_faces)
{
  for (int i=0;i<3;i++){
    // prepare field name
    std::string field_name;
    int cx, cy, cz, delta;
    if (face_center) {
      cx = 0; cy = 0; cz = 0;
      delta = (exterior_faces) ? 1 : -1;
    } else {
      if (exterior_faces){
	cx = 1; cy = 1; cz = 1;
	delta = -1;
      } else {
	cx = -1; cy = -1; cz = -1;
	delta = 1;
      }
    }

    if (i == 0){
      cx += delta;
      field_name = field_prefix + "x";
    } else if (i == 1) {
      cy += delta;
      field_name = field_prefix + "y";
    } else {
      cz += delta;
      field_name = field_prefix + "z";
    }

    // reserve/allocate field
    EnzoTempFieldUtils::prep_reused_temp_field(field, field_name, cx, cy, cz);

    grouping.add(field_name,group_name);
  }
}

//----------------------------------------------------------------------

EnzoConstrainedTransport::EnzoConstrainedTransport(Block *block,
						   int num_partial_timesteps)
  : block_(block),
    num_partial_timesteps_(num_partial_timesteps),
    partial_timestep_index_(0)
{
  ASSERT("EnzoConstrainedTransport", "num_partial_timesteps must be positive",
	 num_partial_timesteps_ > 0);

  if (num_partial_timesteps_ != 2){
    ERROR("EnzoConstrainedTransport",
          "This machinery hasn't been tested for cases when "
          "num_partial_timesteps!=2.");
  }

  // setup the group of permanent magnetic fields
  bfieldi_group_.add("bfieldi_x", "bfield");
  bfieldi_group_.add("bfieldi_y", "bfield");
  bfieldi_group_.add("bfieldi_z", "bfield");

  Field field = block->data()->field();

  // reserve/allocate fields for weight fields
  // these are face-centered fields that store the upwind/downwind direction
  prep_temp_vector_grouping_(field, "weight", weight_group_,
			     "temp_weight_", true, false);

  // allocate temporary efield fields
  // reserve/allocate fields for edge-centered electric fields
  prep_temp_vector_grouping_(field, "efield", efield_group_,
			     "temp_efield_",false,false);

  // reserve/allocate cell-centered e-field
  center_efield_name_ = "center_efield";
  EnzoTempFieldUtils::prep_reused_temp_field(field, center_efield_name_,
					     0, 0, 0);

  if (num_partial_timesteps_ > 1){
    // reserve allocate temporary interface bfields fields (includes the
    // exterior faces of the grid)
    prep_temp_vector_grouping_(field, "bfield", temp_bfieldi_group_,
			       "temp_bfieldi_",true, true);
  }
}

//----------------------------------------------------------------------

EnzoConstrainedTransport::~EnzoConstrainedTransport()
{
  using EnzoTempFieldUtils::deallocate_grouping_fields;

  Field field = block_->data()->field();
  // deallocate electric fields
  std::vector<std::string> efield_group_names{"efield"};
  deallocate_grouping_fields(field, efield_group_names, efield_group_);
  field.deallocate_temporary(field.field_id(center_efield_name_));

  // deallocate weights
  std::vector<std::string> weight_group_names{"weight"};
  deallocate_grouping_fields(field, weight_group_names, weight_group_);

  if (num_partial_timesteps_ > 1){
    // deallocate the temporary longitudinal bfields
    std::vector<std::string> bfieldi_group_names{"bfield"};
    deallocate_grouping_fields(field, bfieldi_group_names,
			       temp_bfieldi_group_);
  }
}

//----------------------------------------------------------------------

void EnzoConstrainedTransport::update_refresh(Refresh* refresh)
{
  FieldDescr * field_descr = cello::field_descr();
  ASSERT("EnzoConstrainedTransport::update_refresh",
	 ("There must be face-centered permanent fields called \"bfieldi_x\","
	  "\"bfield_y\" and \"bfield_z\"."),
	 (field_descr->is_field("bfieldi_x") &&
	  field_descr->is_field("bfieldi_y") &&
	  field_descr->is_field("bfieldi_z")));

  std::vector<std::string> names = {"bfieldi_x", "bfieldi_y", "bfieldi_z"};
  std::vector<std::string> axes = {"x", "y", "z"};
  for (std::size_t i = 0; i < 3; i++){
    std::string name = names[i];
    // first check that field exists
    ASSERT1("EnzoConstrainedTransport::update_refresh",
	    "There must be face-centered permanent fields called \"%s\"",
	    name.c_str(), field_descr->is_field(name));

    int field_id = field_descr->field_id(name);
    // next check the centering of the field
    int centering[3] = {0, 0, 0};
    field_descr->centering(field_id, &centering[0], &centering[1],
			   &centering[2]);
    for (std::size_t j = 0; j<3; j++){
      if (j!=i){
	ASSERT2("EnzoConstrainedTransport::update_refresh",
		"The \"%s\" field must be cell-centered along the %s-axis",
		name.c_str(), axes[j].c_str(), centering[j] == 0);
      } else {
	ASSERT2("EnzoConstrainedTransport::update_refresh",
		"The \"%s\" field must be face-centered along the %s-axis",
		name.c_str(), axes[j].c_str(), centering[j] == 1);
      }
    }

    // finally add the field to refresh
    refresh->add_field(field_id);
  }
}

//----------------------------------------------------------------------

void EnzoConstrainedTransport::increment_partial_timestep() throw()
{
  partial_timestep_index_++;
  ASSERT2("EnzoConstrainedTransport::increment_partial_timestep",
	  ("This should not be called more than %d times to transition "
	   "between the %d partial timestep(s)"),
	  num_partial_timesteps_-1, num_partial_timesteps_,
	  num_partial_timesteps_ > partial_timestep_index_);
}

//----------------------------------------------------------------------

void EnzoConstrainedTransport::correct_reconstructed_bfield
(Grouping &l_group, Grouping &r_group, int dim, int stale_depth)
{
  Grouping *cur_bfieldi_group;
  if (partial_timestep_index_== 0){
    cur_bfieldi_group = &bfieldi_group_;
  } else if (partial_timestep_index_== 1){
    cur_bfieldi_group = &temp_bfieldi_group_;
  } else {
    ERROR1("EnzoConstrainedTransport::correct_reconstructed_bfield",
	   "Unsure how to handle current partial_timestep_index_ of %d",
	   partial_timestep_index_);
  }

  EnzoFieldArrayFactory array_factory(block_, stale_depth);
  EFlt3DArray bfield, l_bfield,r_bfield;
  bfield = array_factory.interior_bfieldi(*cur_bfieldi_group, dim);
  l_bfield = array_factory.reconstructed_field(l_group, "bfield", dim, dim);
  r_bfield = array_factory.reconstructed_field(r_group,"bfield", dim, dim);

  // All 3 array objects are the same shape
  for (int iz=0; iz<bfield.shape(0); iz++) {
    for (int iy=0; iy<bfield.shape(1); iy++) {
      for (int ix=0; ix<bfield.shape(2); ix++) {
	l_bfield(iz,iy,ix) = bfield(iz,iy,ix);
	r_bfield(iz,iy,ix) = bfield(iz,iy,ix);
      }
    }
  }
  // Equivalently: l_bfield.subarray() = bfield;
  //               r_bfield.subarray() = bfield;
}

//----------------------------------------------------------------------

void EnzoConstrainedTransport::identify_upwind(Grouping &flux_group, int dim,
					       int stale_depth)
{
  EnzoFieldArrayFactory array_factory(block_, stale_depth);

  //  - Currently, weight is set to 1.0 if upwind is in positive direction of
  //    the current dimension, 0 if upwind is in the negative direction of the
  //    current dimension, or 0.5 if there is no upwind direction
  //  - At present, the weights are unnecessary (the same information is
  //    encoded in density flux to figure out this information). However, this
  //    functionallity is implemented in case we decide to adopt the weighting
  //    scheme from Athena++, which requires knowledge of the reconstructed
  //    densities.

  EFlt3DArray density_flux, weight_field;
  density_flux = array_factory.from_grouping(flux_group, "density", 0);
  weight_field = array_factory.from_grouping(weight_group_, "weight", dim);

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

//----------------------------------------------------------------------

void EnzoConstrainedTransport::update_all_bfield_components
(Grouping &cur_prim_group, Grouping &xflux_group, Grouping &yflux_group,
 Grouping &zflux_group, Grouping &out_centered_bfield_group, enzo_float dt,
 int stale_depth)
{
  Grouping *out_bfieldi_group;
  if ((partial_timestep_index_== 1) || (num_partial_timesteps_ == 1)){
    out_bfieldi_group = &bfieldi_group_;
  } else if (partial_timestep_index_== 0){
    out_bfieldi_group = &temp_bfieldi_group_;
  } else {
    ERROR1("EnzoConstrainedTransport::correct_reconstructed_bfield",
	   "Unsure how to handle current partial_timestep_index_ of %d",
	   partial_timestep_index_);
  }

  // First, compute the edge-centered Electric fields (each time, it uses
  // the current integrable quantities
  EnzoConstrainedTransport::compute_all_edge_efields
    (block_, cur_prim_group, xflux_group, yflux_group, zflux_group,
     center_efield_name_, efield_group_, weight_group_, stale_depth);

  // Update longitudinal B-field (add source terms of constrained transport)
  for (int dim = 0; dim<3; dim++){
    EnzoConstrainedTransport::update_bfield
      (block_, dim, efield_group_, bfieldi_group_, *out_bfieldi_group,
       dt, stale_depth);
  }

  // Finally, update cell-centered B-field
  for (int dim = 0; dim<3; dim++){
    EnzoConstrainedTransport::compute_center_bfield
      (block_, dim, out_centered_bfield_group, *out_bfieldi_group,
       stale_depth);
  }
}

//----------------------------------------------------------------------

void EnzoConstrainedTransport::compute_center_efield
(Block *block, int dim, std::string center_efield_name, Grouping &prim_group,
 int stale_depth)
{
  // Load the E-field
  EnzoFieldArrayFactory array_factory(block,stale_depth);
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

  for (int iz=0; iz<efield.shape(0); iz++) {
    for (int iy=0; iy<efield.shape(1); iy++) {
      for (int ix=0; ix<efield.shape(2); ix++) {
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
void EnzoConstrainedTransport::compute_edge_efield
(Block *block, int dim, std::string center_efield_name, Grouping &efield_group,
 Grouping &jflux_group, Grouping &kflux_group, Grouping &weight_group,
 int stale_depth)
{

  EnzoPermutedCoordinates coord(dim);
  // determine components of j and k unit vectors:
  int j_x, j_y, j_z, k_x, k_y, k_z;
  coord.j_unit_vector(j_z, j_y, j_x);
  coord.k_unit_vector(k_z, k_y, k_x);

  // Initialize Cell-Centered E-fields
  EFlt3DArray Ec, Ec_jp1, Ec_kp1, Ec_jkp1;
  EnzoFieldArrayFactory array_factory(block,stale_depth);
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

void EnzoConstrainedTransport::compute_all_edge_efields
  (Block *block, Grouping &prim_group, Grouping &xflux_group,
   Grouping &yflux_group, Grouping &zflux_group, std::string center_efield_name,
   Grouping &efield_group, Grouping &weight_group, int stale_depth)
{

  for (int i = 0; i < 3; i++){
    EnzoConstrainedTransport::compute_center_efield
      (block, i, center_efield_name, prim_group, stale_depth);

    Grouping *jflux_group;
    Grouping *kflux_group;
    if (i == 0){
      jflux_group = &yflux_group;
      kflux_group = &zflux_group;
    } else if (i==1){
      jflux_group = &zflux_group;
      kflux_group = &xflux_group;
    } else {
      jflux_group = &xflux_group;
      kflux_group = &yflux_group;
    }

    EnzoConstrainedTransport::compute_edge_efield
      (block, i, center_efield_name, efield_group, *jflux_group, *kflux_group,
       weight_group, stale_depth);
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
void EnzoConstrainedTransport::update_bfield(Block *block, int dim,
					     Grouping &efield_group,
					     Grouping &cur_bfieldi_group,
					     Grouping &out_bfieldi_group,
					     enzo_float dt, int stale_depth)
{
  EnzoPermutedCoordinates coord(dim);

  // compute the ratios of dt to the widths of cells along j and k directions
  EnzoBlock *enzo_block = enzo::block(block);
  enzo_float dtdj = dt/enzo_block->CellWidth[coord.j_axis()];
  enzo_float dtdk = dt/enzo_block->CellWidth[coord.k_axis()];

  // The following comments all assume that we are talking about unstaled
  // region (and that we have dropped all staled cells)
  EnzoFieldArrayFactory array_factory(block,stale_depth);
  
  // Load edge centered efields 
  EFlt3DArray E_j, ej_Lk, ej_Rk, E_k, ek_Lj, ek_Rj;
  E_j = array_factory.from_grouping(efield_group, "efield", coord.j_axis());
  E_k = array_factory.from_grouping(efield_group, "efield", coord.k_axis());

  // Load interface bfields field (includes exterior faces)
  EFlt3DArray cur_bfield, bcur, out_bfield, bout;
  cur_bfield = array_factory.from_grouping(cur_bfieldi_group, "bfield", dim);
  out_bfield = array_factory.from_grouping(out_bfieldi_group, "bfield", dim);

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
void EnzoConstrainedTransport::compute_center_bfield(Block *block, int dim,
						     Grouping &bfieldc_group,
						     Grouping &bfieldi_group,
						     int stale_depth)
{
  EnzoPermutedCoordinates coord(dim);
  EnzoFieldArrayFactory array_factory(block,stale_depth);

  // Load cell-centerd field
  EFlt3DArray b_center = array_factory.from_grouping(bfieldc_group, "bfield",
						     coord.i_axis());
  // Load Face-centered fields
  EFlt3DArray bi_left = array_factory.from_grouping(bfieldi_group, "bfield",
						    coord.i_axis());
  // Get the view of the Face-center field that starting from i=1
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
