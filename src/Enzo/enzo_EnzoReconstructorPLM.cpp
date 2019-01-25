#include "cello.hpp"
#include "enzo.hpp"

inline enzo_float sign(enzo_float val){
  // implements sign function https://stackoverflow.com/questions/1903954
  return (0.0 < val) - (val < 0.0);
}

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

inline enzo_float monotized_difference(enzo_float vm1, enzo_float v,
				       enzo_float vp1)
{
  // Adapted from Enzo's Rec_PLM.C
  // Original Enzo code applied a limiter
  enzo_float dv_l = (v-vm1);
  enzo_float dv_r = (vp1-v);
  enzo_float dv_c = 0.25*(vp1-vm1);

  // Apply monoticity constraint
  return sign(dv_c)*Min(fabs(dv_l),fabs(dv_r),fabs(dv_c));
}

//----------------------------------------------------------------------

// helper function that sets both:
//   - left reconstructed value at the interface between first and second
//     cells along dimension i, and
//   - right reconstructed value at the interface between the second-to-last
//     and last cells along dimension i
// to zero. Other reconstruction methods (e.g. Nearest Neighbor) may actually
// be to reconstruct values at these locations.
//
// Arguments:
//  - wl and wr are both face-centered along dimension i. These arrays do not
//    include faces on the exterior of the grid.
//  - di{dim}, where {dim} is x,y,z indicates the amount by which i{dim}
//    changes if we move in the i direction. For example, if the ith
//    dimension aligns with z: dix=0, diy=0, diz = 1
void zero_edge_values_(EnzoArray<enzo_float> &wl, EnzoArray<enzo_float> &wr,
		       int diz, int diy, int dix, enzo_float prim_floor)
{
  // Number of cells in the grid along each dimension
  int mz = wl.length_dim2() + diz;
  int my = wl.length_dim1() + diy;
  int mx = wl.length_dim0() + dix;

  // mk,mj equal number of grid cells along k,j directions 
  // In both cases, iterate: k = 0  to (but not including) k = mk
  //                         j = 0  to (but not including) j = mj 
  // For wl, iterate:        i = 0  to (but not including) i = 1
  for (int iz=0; iz<(mz*(1-diz) + diz); iz++) {
    for (int iy=0; iy<(my*(1-diy) + diy); iy++) {
      for (int ix=0; ix<(mx*(1-dix) + dix); ix++) {
        wl(iz,iy,ix) = prim_floor;
      }
    }
  }

  // For wr, iterate:        i = (mi-2) up to (but not including) i = (mi-1)
  for (int iz=(diz*(mz-2)); iz<(mz-diz); iz++) {
    for (int iy=(diy*(my-2)); iy<(my-diy); iy++) {
      for (int ix=(dix*(mx-2)); ix<(mx-dix); ix++) {
        wr(iz,iy,ix) = prim_floor;
      }
    }
  }


  // The above may is missing values
  // Once it is fixed, the following can be removed.
  for (int iz=0; iz<wl.length_dim2(); iz++) {
    for (int iy=0; iy<wl.length_dim1(); iy++) {
      for (int ix=0; ix<wl.length_dim0(); ix++) {
        wl(iz,iy,ix) = prim_floor;
  	wr(iz,iy,ix) = prim_floor;
      }
    }
  }
  
}

//----------------------------------------------------------------------
void EnzoReconstructorPLM::reconstruct_interface (Block *block,
						  Grouping &prim_group,
						  Grouping &priml_group,
						  Grouping &primr_group,
						  int dim,
						  EnzoEquationOfState *eos)
{
  std::vector<std::string> prim_group_names = EnzoMethodVlct::prim_group_names;
  EnzoFieldArrayFactory array_factory;

  // compute which index needs to be changed to advance in direction of dim
  int dix = (dim == 0) ? 1 : 0;
  int diy = (dim == 1) ? 1 : 0;
  int diz = (dix == diy) ? 1 : 0;

  // unecessary values are computed for the inside faces of outermost ghost zone
  for (unsigned int group_ind=0;group_ind<prim_group_names.size(); group_ind++){

    // load group name and number of fields in the group
    std::string group_name = prim_group_names[group_ind];
    int num_fields = prim_group.size(group_name);

    // Handle possibility of having a density/pressure floor
    enzo_float prim_floor =0;
    bool use_floor = false;
    if (group_name == "density"){
      prim_floor = eos->get_density_floor();
      use_floor=true;
    } else if (group_name == "pressure"){
      prim_floor = eos->get_pressure_floor();
      use_floor=true;
    }

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){
      // Cast the problem as reconstructing values at:
      //   wl(k, j, i+3/2) and wr(k,j,i+1/2)

      // Prepare cell-centered arrays
      // define:   wc_left(k,j,i)   -> w(k,j,i)
      //           wc_center(k,j,i) -> w(k,j,i+1)
      //           wc_right(k,j,i)  -> w(k,j,i+2)
      EnzoArray<enzo_float> wc_left, wc_center, wc_right;
      array_factory.load_grouping_field(block, prim_group, group_name,
					field_ind, wc_left);
      wc_center.initialize_subarray(wc_left, diz, wc_left.length_dim2(),
				    diy, wc_left.length_dim1(),
				    dix, wc_left.length_dim0());
      wc_right.initialize_subarray(wc_left, 2*diz, wc_left.length_dim2(),
				   2*diy, wc_left.length_dim1(),
				   2*dix, wc_left.length_dim0());

      
      // Prepare face-centered arrays
      // define:   wl_offset(k,j,i)-> wl(k,j,i+3/2)
      //           wr(k,j,i)       -> wr(k,j,i+1/2)
      EnzoArray<enzo_float> wr, wl, wl_offset;
      array_factory.load_temp_interface_grouping_field(block, primr_group,
						       group_name, field_ind,
						       wr, dim != 0, dim != 1,
						       dim != 2);
      array_factory.load_temp_interface_grouping_field(block, priml_group,
						       group_name, field_ind,
						       wl, dim != 0, dim != 1,
						       dim != 2);
      wl_offset.initialize_subarray(wl, diz, wl.length_dim2(),
				    diy, wl.length_dim1(),
				    dix, wl.length_dim0());

      // At the interfaces between the first and second cell (second-to-
      // last and last cell), along a given axis, set the reconstructed
      // left (right) interface value to prim_floor (or 0)
      // (the use of wl instead of wl_offset is intentional)
      zero_edge_values_(wl, wr, diz, diy, dix, prim_floor);
      
      // Iteration Limits are compatible with a 2D and 3D grid
      // Need to reconstruct and compute fluxes at k=j=0 and k=kmax-1, j=jmax-1
      // so that we can perform constrained transport
      // Iterate over: k = 0, 1, ..., kmax - 1 
      //               j = 0, 1, ..., jmax - 1
      //               i = 0, 1, ..., imax - 3

      for (int iz=0; iz<wc_left.length_dim2()-2*diz; iz++) {
	for (int iy=0; iy<wc_left.length_dim1()-2*diy; iy++) {
	  for (int ix=0; ix<wc_left.length_dim0()-2*dix; ix++) {

	    // compute monotized difference
	    enzo_float val = wc_center(iz,iy,ix);
	    enzo_float dv = monotized_difference(wc_left(iz,iy,ix), val,
						 wc_right(iz,iy,ix));
	    enzo_float left_val, right_val;
	    right_val = use_floor ? std::max(val-dv,prim_floor) : val-dv;
	    left_val = use_floor ? std::max(val+dv,prim_floor) : val+dv;

	    // face centered fields: index i corresponds to the value at i-1/2
	    wr(iz,iy,ix) = left_val;
	    wl_offset(iz,iy,ix) = right_val;
	  }
	}
      }
    }
  }
}
