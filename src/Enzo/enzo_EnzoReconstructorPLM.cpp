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

// Need to place Floor Values!

void EnzoReconstructorPLM::reconstruct_interface (Block *block,
						  Grouping &prim_group,
						  Grouping &priml_group,
						  Grouping &primr_ids, int dim)
{
  Field field = block->data()->field();

  // Retrieve iteration limits
  const int id = prim_ids[0];

  // cell-centered iteration dimensions
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);

  // compute face-centered dimensions and offset between neigboring elements
  // along dim (offset is the same for cell-centered and face-centered arrays)
  int fc_mx = mx;
  int fc_my = my;
  int fc_mz = mz;

  int offset;
  if (dim == 0) {
    fc_mx++;
    offset = 1;
  } else if (dim == 1) {
    fc_my++;
    offset = mx;
  } else {
    fc_mz++;
    offset = mx+my;
  }

  // number of fields for which we will perform reconstruction
  int nfields = prim_ids.size();

  // In the current implementation that follows, unecessary values are computed
  // for face just interior to the outermost ghost zone
  for (int group_ind=0; group_ind<num_cons_names; group_ind++){

    // load group name and number of fields in the group
    std::string group_name = prim_group_names[group_ind];
    int num_fields = prim_group.size(group_name);

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      // Load in the primitives
      enzo_float *w = load_grouping_field_(field, &prim_group, group_name,
					   field_ind);
      enzo_float *wl = load_grouping_field_(field, &priml_group, group_name,
					    field_ind);
      enzo_float *wr = load_grouping_field_(field, &primr_group, group_name,
					   field_ind);

      for (int iz=0; iz<fc_mz; iz++) {
	for (int iy=0; iy<fc_my; iy++) {
	  for (int ix=0; ix<fc_mx; ix++) {

	    // compute the index of the cell-centered and face-centered
	    int fc_i = ix + fc_mx*(iy + fc_my*iz);

	    if ((ix == 0) || (iy == 0) || (iz == 0)){
	      // At the interface between the first and second cell along a
	      // given axis, set the reconstructed left interface value to 0
	      wl[fc_i+offset] = 0;
	      continue;
	    } else if ((ix == fc_mx-1) || (ix == fc_mx-1) || (ix == fc_mx-1)){
	      // At the interface between the second-to-last and last cell
	      // along a given axis, set the reconstructed right interface
	      // value to 0
	      wr[fc_i] = 0;
	      continue;
	    }
	    
	    // compute the index of the cell-centered quantity
	    int cc_i = ix + mx*(iy + my*iz);

	    enzo_float val = w[cc_i];

	    // compute monotized difference
	    enzo_float dv = monotized_difference(w[cc_i-offset], val,
						 w[cc_i+offset]);

	    // for face centered values, index i corresponds to the value at
	    // i-1/2
	    wr[fc_i] = val-dv;
	    wl[fc_i+offset] = val+dv;
	  
	  }
	}
      }
    }
  }
  
}
