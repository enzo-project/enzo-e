// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructorPLM.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implements the EnzoReconstructorPLM class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

inline enzo_float sign(enzo_float val){
  // implements sign function https://stackoverflow.com/questions/1903954
  return (0.0 < val) - (val < 0.0);
}

//----------------------------------------------------------------------

inline enzo_float Min(enzo_float a, enzo_float b, enzo_float c)
{
  // taken from Enzo's ReconstructionRoutines.h
  if (a<b) {
    if (c<a)
      return c;
    else 
      return a;
  } else {
    if (c<b)
      return c;
    else 
      return b;
  }
}

//----------------------------------------------------------------------

inline enzo_float monotized_difference(enzo_float vm1, enzo_float v,
				       enzo_float vp1)
{
  // The limiter is expected to return 0.5 times the limted slope

  // Stone & Gardiner (2008) misleadingly mention a minmod limiter while
  // discussing this limiter but the textbook they cite (LeVeque 2002) refers
  // to it as a monotized central-difference limiter). The limiter shown in the
  // paper differs from the one used in Athena for primitives

  // Always return 0 if sign(dv_l) != sign(dv_r)
  // When sign(dv_l) == sign(dv_r), there is a difference in the limter used by
  // Enzo (which matches the description of LeVeque (2002)) and the limiter
  // used in the c version of Athena (for reconstructing primitives).
  //   - Enzo would return:
  //       sign(dv_c)*Min(fabs(dv_l),fabs(dv_r),fabs(dv_c));
  //   - Athena would compute dv_g = dv_l*dv_r/(dv_l+dv_r) and would return:
  //       sign(dv_c)*Min(fabs(dv_l),fabs(dv_r),fabs(dv_c),fabs(dv_g))
  //     If this is used and if a resulting reconstructed value exceeds the
  //     cell-centered value of one of its direct neighbors, that reconstructed
  //     value MUST be rounded to the closest neighboring value

  // Current implementation adapted from Enzo's Rec_PLM.C
  enzo_float dv_l = (v-vm1);
  enzo_float dv_r = (vp1-v);
  enzo_float dv_c = 0.25*(vp1-vm1);

  return 0.5*(sign(dv_l)+sign(dv_r))*Min(fabs(dv_l),fabs(dv_r),fabs(dv_c));
}

//----------------------------------------------------------------------

/// helper function that sets both left reconstructed value at the interface
/// between first and second cells along dimension i, and right reconstructed
/// value at the interface between the second-to-last and last cells along
/// dimension i to zero (or floor). Other reconstruction methods (e.g. Nearest
/// Neighbor) may actually be to reconstruct values at these locations.
///
/// @param wl,wr face-centered arrays along dimension dim. These arrays do not
///  include faces on the exterior of the grid.
/// @param prim_floor value to set the edges equal to
/// @dim the dimension transverse normal to the edge (the dimension along which
///  we reconstruct values
void zero_edge_values_(EFlt3DArray &wl, EFlt3DArray &wr, enzo_float prim_floor,
		       int dim)
{
  // Implementation would be simplified if we didn't have to specify the stop
  // index of the slice 
  
  if (dim == 0){
    // reconstructed along x-axis
    wl.subarray(CSlice( 0,wl.shape(0)),
		CSlice( 0,wl.shape(1)),
		CSlice( 0,          1)) = prim_floor;
    wr.subarray(CSlice( 0,wr.shape(0)),
		CSlice( 0,wr.shape(1)),
		CSlice(-1,wr.shape(2))) = prim_floor;
  } else if (dim == 1){
    // reconstructed along y-axis
    wl.subarray(CSlice( 0,wl.shape(0)),
		CSlice( 0,          1),
		CSlice( 0,wl.shape(2))) = prim_floor;
    wr.subarray(CSlice( 0,wr.shape(0)),
		CSlice(-1,wr.shape(1)),
		CSlice( 0,wr.shape(2))) = prim_floor;
  } else {
    // reconstructed along y-axis
    wl.subarray(CSlice( 0,          1),
		CSlice( 0,wl.shape(1)),
		CSlice( 0,wl.shape(2))) = prim_floor;
    wr.subarray(CSlice(-1,wr.shape(0)),
		CSlice( 0,wr.shape(1)),
		CSlice( 0,wr.shape(2))) = prim_floor;
  }
}

//----------------------------------------------------------------------

void EnzoReconstructorPLM::reconstruct_interface (Block *block,
						  Grouping &prim_group,
						  Grouping &priml_group,
						  Grouping &primr_group,
						  int dim,
						  EnzoEquationOfState *eos,
						  int stale_depth)
{
  std::vector<std::string> group_names = this->group_names_;

  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EnzoPermutedCoordinates coord(dim);

  // unecessary values are computed for the inside faces of outermost ghost zone
  for (unsigned int group_ind=0; group_ind<group_names.size(); group_ind++){

    // load group name and number of fields in the group
    std::string group_name = group_names[group_ind];
    int num_fields = prim_group.size(group_name);

    // Handle possibility of having a density floor
    enzo_float prim_floor =0;
    bool use_floor = false;
    if (group_name == "density"){
      prim_floor = eos->get_density_floor();
      use_floor=true;
    } else if (group_name == "internal_energy"){
      // This is just used to initialize the external values of the grid, to a
      // known allowed internal energy. (This may not be necessary)
      // The floor on the internal_energy is checked separately at the end of
      // the function
      prim_floor = (eos->get_pressure_floor() /
		    (eos->get_density_floor() * (eos->get_gamma()-1.)));
      use_floor = false;
    }

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){
      // Cast the problem as reconstructing values at:
      //   wl(k, j, i+3/2) and wr(k,j,i+1/2)

      // Prepare cell-centered arrays
      // define:   wc_left(k,j,i)   -> w(k,j,i)
      //           wc_center(k,j,i) -> w(k,j,i+1)
      //           wc_right(k,j,i)  -> w(k,j,i+2)
      EFlt3DArray wc_left, wc_center, wc_right;
      wc_left = array_factory.from_grouping(prim_group, group_name, field_ind);
      wc_center = coord.left_edge_offset(wc_left, 0, 0, 1);
      wc_right  = coord.left_edge_offset(wc_left, 0, 0, 2);

      // Prepare face-centered arrays
      // define:   wl_offset(k,j,i)-> wl(k,j,i+3/2)
      //           wr(k,j,i)       -> wr(k,j,i+1/2)
      EFlt3DArray wr, wl, wl_offset;
      wr = array_factory.reconstructed_field(primr_group, group_name, field_ind,
					     dim);
      wl = array_factory.reconstructed_field(priml_group, group_name, field_ind,
					     dim);
      wl_offset = coord.left_edge_offset(wl, 0, 0, 1);


      // Now that we have introduced immediate_stale depth the following should
      // not be necessary
      // At the interfaces between the first and second cell (second-to-
      // last and last cell), along a given axis, set the reconstructed
      // left (right) interface value to prim_floor (or 0)
      // (the use of wl instead of wl_offset is intentional)
      zero_edge_values_(wl, wr, prim_floor, dim);
      
      // Iteration Limits are compatible with a 2D and 3D grid
      // Need to reconstruct and compute fluxes at k=j=0 and k=kmax-1, j=jmax-1
      // so that we can perform constrained transport
      // Iterate over: k = 0, 1, ..., kmax - 1 
      //               j = 0, 1, ..., jmax - 1
      //               i = 0, 1, ..., imax - 3

      for (int iz=0; iz<wc_right.shape(0); iz++) {
	for (int iy=0; iy<wc_right.shape(1); iy++) {
	  for (int ix=0; ix<wc_right.shape(2); ix++) {

	    // compute monotized difference
	    enzo_float val = wc_center(iz,iy,ix);
	    enzo_float dv = monotized_difference(wc_left(iz,iy,ix), val,
						 wc_right(iz,iy,ix));
	    enzo_float left_val, right_val;
	    // It's worth mentioning that in the original Athena c code, rather
	    // than placing floors on a few select quantites, they instead
	    // checked that all the reconstructed quantities lied between the
	    // left and right cell-centered quantities

	    if (use_floor) {
	      right_val = EnzoEquationOfState::apply_floor(val-dv,prim_floor);
	      left_val = EnzoEquationOfState::apply_floor(val+dv,prim_floor);
	    } else {
	      right_val = val-dv;
	      left_val = val+dv;
	    }

	    // If using the Athena limiter, need to also use:
	    //right_val = std::max(std::min(val,wc_left(iz,iy,ix)),val-dv);
	    //right_val = std::min(std::max(val,wc_left(iz,iy,ix)),right_val);
	    //left_val = std::max(std::min(val,wc_right(iz,iy,ix)),val+dv);
	    //left_val = std::min(std::max(val,wc_right(iz,iy,ix)),left_val);

	    // face centered fields: index i corresponds to the value at i-1/2
	    wr(iz,iy,ix) = right_val;
	    wl_offset(iz,iy,ix) = left_val;
	  }
	}
      }
    }
  }


  // apply the floor to the internal energy (increment stale_depth by 1 since
  // this has an immediate stale depth of 1)
  eos->apply_floor_to_internal_energy(block, priml_group, stale_depth+1);
  eos->apply_floor_to_internal_energy(block, primr_group, stale_depth+1);
}
