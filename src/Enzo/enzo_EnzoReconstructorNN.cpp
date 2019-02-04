#include "cello.hpp"
#include "enzo.hpp"


void EnzoReconstructorNN::reconstruct_interface (Block *block,
						  Grouping &prim_group,
						  Grouping &priml_group,
						  Grouping &primr_group,
						  int dim,
						  EnzoEquationOfState *eos)
{
  std::vector<std::string> prim_group_names = EnzoMethodVlct::prim_group_names;

  EnzoFieldArrayFactory array_factory(block);

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
    // cell-centered values should already satisfy the floor

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      // define wc_offset(k,j,i) -> wc(k,j,i+1)
      EFlt3DArray wc = array_factory.from_grouping(prim_group, group_name,
						    field_ind);
      EFlt3DArray wc_offset = wc.subarray(diz, wc.length_dim2(),
					  diy, wc.length_dim1(),
					  dix, wc.length_dim0());

      EFlt3DArray wr, wl;
      // face centered: (k,j,i) -> (k,j,i+1/2)
      wr = array_factory.load_temp_interface_grouping_field(primr_group,
							    group_name,
							    field_ind, dim != 0,
							    dim != 1, dim != 2);
      wl = array_factory.load_temp_interface_grouping_field(priml_group,
							    group_name,
							    field_ind, dim != 0,
							    dim != 1, dim != 2);

      // Could be more efficient
      for (int iz=0; iz<wc.length_dim2()-diz; iz++) {
	for (int iy=0; iy<wc.length_dim1()-diy; iy++) {
	  for (int ix=0; ix<wc.length_dim0()-dix; ix++) {
	      wl(iz,iy,ix) = wc(iz,iy,ix);
	      wr(iz,iy,ix) = wc_offset(iz,iy,ix);
	  }
	}
      }
    }
  }
}
