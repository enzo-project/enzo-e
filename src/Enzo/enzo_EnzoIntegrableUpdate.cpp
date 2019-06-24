// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoIntegrableUpdate.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon June 24 2019
/// @brief    [\ref Enzo] Implementation of EnzoIntegrableUpdate.

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::setup_lut_(){
  EnzoCenteredFieldRegistry registry;
  lut_ = registry.prepare_advection_lut(integrable_groups_,
					conserved_start_, conserved_stop_,
					specific_start_, specific_stop_,
				        other_start_, other_stop_, nfields_,
					flagged_groups_);

  ASSERT("EnzoMethodMHDVlct::update_quantities_",
	 "not equipped to handle fields classified as other",
	 other_start_ == other_stop_);
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::update_quantities
  (Block *block, Grouping &initial_integrable_group,
   Grouping &xflux_group, Grouping &yflux_group, Grouping &zflux_group,
   Grouping &out_integrable_group, Grouping &out_conserved_passive_scalar,
   EnzoEquationOfState *eos, double dt, int stale_depth)
{

  // Update passive scalars, it doesn't currently support renormalizing to 1
  update_passive_scalars_(block, initial_integrable_group,
			  xflux_group, yflux_group, zflux_group,
			  out_conserved_passive_scalar,  dt, stale_depth);

  EnzoCenteredFieldRegistry registry;
  EFlt3DArray *cur_prim_arrays, *out_prim_arrays;
  cur_prim_arrays = registry.load_array_of_fields(block, lut_, nfields_,
						  initial_integrable_group,
						  stale_depth);
  out_prim_arrays = registry.load_array_of_fields(block, lut_, nfields_,
						  out_integrable_group,
						  stale_depth);
  
  EFlt3DArray *xflux_arrays, *yflux_arrays, *zflux_arrays;
  // Note: No need to identify x-flux, y-flux, z-flux as reconstructed
  //       since the shape of the arrays are appropriately registerred
  //       (Additionally, if we specified dim, then the velocities saved to
  //        axes i, j, k would be different for each component)
  xflux_arrays = registry.load_array_of_fields(block, lut_, nfields_,
					       xflux_group, stale_depth);
  yflux_arrays = registry.load_array_of_fields(block, lut_, nfields_,
					       yflux_group, stale_depth);
  zflux_arrays = registry.load_array_of_fields(block, lut_, nfields_,
					       zflux_group, stale_depth);

  enzo_float *cur_prim, *dU;
  cur_prim = new enzo_float[nfields_];
  dU = new enzo_float[nfields_];
  
  // For now, not having density floor affect momentum or total energy density
  enzo_float density_floor = eos->get_density_floor();

  // cell-centered grid dimensions
  int mz = cur_prim_arrays[lut_.density].shape(0);
  int my = cur_prim_arrays[lut_.density].shape(1);
  int mx = cur_prim_arrays[lut_.density].shape(2);

  EnzoBlock * enzo_block = enzo::block(block);
  enzo_float dtdx = dt/enzo_block->CellWidth[0];
  enzo_float dtdy = dt/enzo_block->CellWidth[1];
  enzo_float dtdz = dt/enzo_block->CellWidth[2];

  for (int iz=1; iz<mz-1; iz++) {
    for (int iy=1; iy<my-1; iy++) {
      for (int ix=1; ix<mx-1; ix++) {

	// load in the fields
	for (int field_ind=0; field_ind<nfields_; field_ind++){
	  cur_prim[field_ind] = cur_prim_arrays[field_ind](iz,iy,ix);
	  dU[field_ind] = (- dtdx * (xflux_arrays[field_ind](iz,iy,ix)
				    - xflux_arrays[field_ind](iz,iy,ix-1))
			   - dtdy * (yflux_arrays[field_ind](iz,iy,ix)
				     - yflux_arrays[field_ind](iz,iy-1,ix))
			   - dtdz * (zflux_arrays[field_ind](iz,iy,ix)
				     - zflux_arrays[field_ind](iz-1,iy,ix)));
	}

	// get the initial density
	enzo_float initial_density = cur_prim[lut_.density];

	// now update the integrable primitives that are conserved
	for (int i = conserved_start_; i < conserved_stop_; i++){
	  out_prim_arrays[i](iz,iy,ix) = cur_prim[i] + dU[i];
	}

	enzo_float new_density = out_prim_arrays[lut_.density](iz,iy,ix);
	// possibly place a floor on new_density.
	new_density = EnzoEquationOfState::apply_floor(new_density,
						       density_floor);
	out_prim_arrays[lut_.density](iz,iy,ix) = new_density;

	// update the specific primitives
	for (int i = specific_start_; i < specific_stop_; i++){
	  out_prim_arrays[i](iz,iy,ix) = (cur_prim[i] * initial_density
					  + dU[i]) /new_density;
	}

	// if relevent, this would be a good place to handle the dual energy
	// formalism ...
      }
    }
  }

  // apply floor to energy
  eos->apply_floor_to_total_energy(block, out_integrable_group, stale_depth+1);

  delete[] cur_prim;         delete[] dU;
  delete[] cur_prim_arrays;  delete[] out_prim_arrays;
  delete[] xflux_arrays;     delete[] yflux_arrays;     delete[] zflux_arrays;
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::update_passive_scalars_
  (Block *block, Grouping &initial_integrable_group, Grouping &xflux_group,
   Grouping &yflux_group, Grouping &zflux_group,
   Grouping &out_conserved_passive_scalar,  double dt, int stale_depth)
{

  // Does not currently handle passively advected scalars that must sum to 1
  // The easiest way to do that would be to divide each group of fields that
  // must sum to 1 into different groupings
  
  EnzoFieldArrayFactory array_factory(block, stale_depth);

  EFlt3DArray cur_density
    = array_factory.from_grouping(initial_integrable_group, "density", 0);
  
  // cell-centered grid dimensions
  int mz = cur_density.shape(0);
  int my = cur_density.shape(1);
  int mx = cur_density.shape(2);

  EnzoBlock * enzo_block = enzo::block(block);
  enzo_float dtdx = dt/enzo_block->CellWidth[0];
  enzo_float dtdy = dt/enzo_block->CellWidth[1];
  enzo_float dtdz = dt/enzo_block->CellWidth[2];

  std::vector<std::string> group_names = this->passive_groups_;

  for (std::size_t group_ind = 0; group_ind<group_names.size(); group_ind++){

    // load group name and number of fields in the group
    std::string group_name = group_names[group_ind];
    int num_fields = initial_integrable_group.size(group_name);

    if (num_fields == 0){
      continue;
    }

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      // load values
      EFlt3DArray cur_specific, out_conserved, xflux, yflux, zflux;
      cur_specific = array_factory.from_grouping(initial_integrable_group,
						 group_name, field_ind);
      out_conserved = array_factory.from_grouping(out_conserved_passive_scalar,
						  group_name, field_ind);
      xflux = array_factory.from_grouping(xflux_group, group_name, field_ind);
      yflux = array_factory.from_grouping(yflux_group, group_name, field_ind);
      zflux = array_factory.from_grouping(zflux_group, group_name, field_ind);

      for (int iz=1; iz<mz-1; iz++) {
	for (int iy=1; iy<my-1; iy++) {
	  for (int ix=1; ix<mx-1; ix++) {

	    out_conserved(iz,iy,ix)
	      = (cur_specific(iz,iy,ix) * cur_density(iz,iy,ix)
		 - dtdx * (xflux(iz,iy,ix) - xflux(iz,iy,ix-1))
		 - dtdy * (yflux(iz,iy,ix) - yflux(iz,iy-1,ix))
		 - dtdz * (zflux(iz,iy,ix) - zflux(iz-1,iy,ix)));

	  }
	}
      }

    }
  }
}
