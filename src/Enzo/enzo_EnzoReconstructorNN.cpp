#include "cello.hpp"
#include "enzo.hpp"


void EnzoReconstructorNN::reconstruct_interface (Block *block,
						  Grouping &prim_group,
						  Grouping &priml_group,
						  Grouping &primr_group,
						  int dim,
						  EnzoEquationOfState *eos)
{
  std::vector<std::string> group_names = this->group_names_;

  EnzoFieldArrayFactory array_factory(block);
  // determine components of i unit vector
  EnzoPermutedCoordinates coord(dim);
  int i_x, i_y, i_z;
  coord.i_unit_vector(i_x,i_y,i_z);

  for (unsigned int group_ind=0; group_ind<group_names.size(); group_ind++){

    // load group name and number of fields in the group
    std::string group_name = group_names[group_ind];
    int num_fields = prim_group.size(group_name);

    // Handle possibility of having a density/pressure floor:
    //  -> No need b/c cell-centered values should already satisfy the floor

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      // define wc_offset(k,j,i) -> wc(k,j,i+1)
      EFlt3DArray wc = array_factory.from_grouping(prim_group, group_name,
						    field_ind);
      EFlt3DArray wc_offset = coord.left_edge_offset(wc, 0, 0, 1);

      EFlt3DArray wr, wl;
      wr = array_factory.reconstructed_field(primr_group, group_name, field_ind,
					     dim);
      wl = array_factory.reconstructed_field(priml_group, group_name, field_ind,
					     dim);

      // Could be more efficient
      for (int iz=0; iz<wc.shape(0)-i_z; iz++) {
	for (int iy=0; iy<wc.shape(1)-i_y; iy++) {
	  for (int ix=0; ix<wc.shape(2)-i_x; ix++) {
	      wl(iz,iy,ix) = wc(iz,iy,ix);
	      wr(iz,iy,ix) = wc_offset(iz,iy,ix);
	  }
	}
      }
    }
  }
}
