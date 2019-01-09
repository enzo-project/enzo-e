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

void EnzoReconstructorPLM::reconstruct_interface (Block *block,
						  Grouping &prim_group,
						  Grouping &priml_group,
						  Grouping &primr_group,
						  int dim,
						  EnzoEquationOfState *eos)
{
  std::vector<std::string> prim_group_names = EnzoMethodVlct::prim_group_names;

  // compute which index needs to be changed to advance in direction of dim
  int dix = (dim == 0) ? 1 : 0;
  int diy = (dim == 1) ? 1 : 0;
  int diz = (dim == 2) ? 1 : 0;

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
      // Load in the primitives
      EnzoArray<enzo_float> w, wl, wr;
      load_grouping_field_(block, prim_group, group_name, field_ind, w);
      load_interface_prim_field_(block, priml_group, group_name, field_ind, wl,
				 dim);
      load_interface_prim_field_(block, priml_group, group_name, field_ind, wr,
				 dim);

      // Get iteration limits
      int zstop = wl.length_dim2();
      int ystop = wl.length_dim1();
      int xstop = wl.length_dim0();

      // if dim points along direction x. Then the following unecessarily sets
      // interface values to 0 when y = 0 || ystop-1. ALSO unecessarily sets
      // interface values to 0 when z = 0 || zmax-1.
      for (int iz=0; iz<zstop; iz++) {
	for (int iy=0; iy<ystop; iy++) {
	  for (int ix=0; ix<xstop; ix++) {

	    // At the interfaces between the first and second cell (second-to-
	    // last and last cell), along a given axis, set the reconstructed
	    // left (right) interface value to 0
	    // (we also set interfaces between cells in the first and last rows
	    //  to zero - this is ok, but unecessary)
	    if ((ix == 0) || (iy == 0) || (iz == 0)){
	      wl(iz+diz,iy+diy,ix+dix) = 0;
	      continue;
	    } else if ((iz == zstop-1) || (iy == ystop-1) || (ix == xstop-1)){
	      wr(iz,iy,ix) = 0;
	      continue;
	    }

	    // compute the index of the cell-centered quantity
	    enzo_float val = w(iz,iy,ix);

	    // compute monotized difference
	    enzo_float dv = monotized_difference(w(iz-diz,iy-diy,ix-dix), val,
						 w(iz+diz,iy+diy,ix+dix));
	    enzo_float left_val, right_val;
	    right_val = use_floor ? val-dv : std::max(val-dv,prim_floor);
	    left_val = use_floor ? val+dv : std::max(val+dv,prim_floor);

	    // face centered fields: index i corresponds to the value at i-1/2
	    wr(iz,iy,ix) = left_val;
	    wl(iz+diz,iy+diy,ix+dix) = right_val;
	  }
	}
      }
    }
  }
}
