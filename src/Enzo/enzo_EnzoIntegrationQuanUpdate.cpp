// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoIntegrationQuanUpdate.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon June 24 2019
/// @brief    [\ref Enzo] Implementation of EnzoIntegrationQuanUpdate.

#include "cello.hpp"
#include "enzo.hpp"
#include <limits>

//----------------------------------------------------------------------

static void append_key_to_vec_
(const str_vec_t &integration_quantities, FieldCat target_cat,
 const bool skip_bfield, std::size_t *density_index, str_vec_t &key_vec)
{
  for (const std::string &name : integration_quantities){
    bool vector_quantity, actively_advected;
    FieldCat category;
    bool success = EnzoCenteredFieldRegistry::quantity_properties
      (name, &vector_quantity, &category, &actively_advected);

    // Sanity Checks:
    ASSERT1("append_key_to_vec_",
	    ("\"%s\" is not registered in EnzoCenteredFieldRegistry"),
	    name.c_str(), success);
    ASSERT1("append_key_to_vec_",
	    ("\"%s\" should not be listed as an integration quantity because "
	     "it is not actively advected."),
	    name.c_str(), actively_advected);

    if (category != target_cat){
      ASSERT1("append_key_to_vec_",
	      ("Can't handle the integration quantity, \"%s\", because it has "
	       "a field category of FieldCat::other"),
	      name.c_str(), category != FieldCat::other);
      continue;
    } else if (skip_bfield && (name == "bfield")){
      continue;
    } else if ((density_index != nullptr) && (name == "density")){
      *density_index = key_vec.size();
    }

    if (vector_quantity){
      key_vec.push_back(name + "_x");
      key_vec.push_back(name + "_y");
      key_vec.push_back(name + "_z");
    } else {
      key_vec.push_back(name);
    }
  }
}

//----------------------------------------------------------------------

EnzoIntegrationQuanUpdate::EnzoIntegrationQuanUpdate
(const str_vec_t& integration_quantities,
 const bool skip_B_update) throw()
{
  // Add conserved quantities to integration_keys_ and identify the index
  // holding the density key
  density_index_ = std::numeric_limits<std::size_t>::max();
  append_key_to_vec_(integration_quantities, FieldCat::conserved, skip_B_update,
                     &density_index_, integration_keys_);
  // Confirm that density is in fact a registered quantity
  ASSERT("EnzoIntegrationQuanUpdate",
	 ("\"density\" must be a registered integration quantity."),
	 density_index_ != std::numeric_limits<std::size_t>::max());
  // Record the first index holding a key for a specific quantity
  first_specific_index_ = integration_keys_.size();
  // Add specific quantities to integration_keys_
  append_key_to_vec_(integration_quantities, FieldCat::specific, skip_B_update,
                     nullptr, integration_keys_);
}

//----------------------------------------------------------------------

void EnzoIntegrationQuanUpdate::clear_dUcons_map
(EnzoEFltArrayMap &dUcons_map, enzo_float value,
 const str_vec_t &passive_list) const noexcept
{
  auto init_arr = [value,&dUcons_map](const std::string& key)
  {
    EFlt3DArray array = dUcons_map.at(key);
    for (int iz=0; iz<array.shape(0); iz++) {
      for (int iy=0; iy<array.shape(1); iy++) {
        for (int ix=0; ix<array.shape(2); ix++) {
          array(iz,iy,ix) = value;
        }
      }
    }
  };

  for (const std::string& key : integration_keys_){ init_arr(key); }
  for (const std::string& key : passive_list){ init_arr(key); }
}

//----------------------------------------------------------------------

void EnzoIntegrationQuanUpdate::accumulate_flux_component
(int dim, double dt, enzo_float cell_width, EnzoEFltArrayMap &flux_map,
 EnzoEFltArrayMap &dUcons_map, int stale_depth,
 const str_vec_t &passive_list) const noexcept
{
  EnzoPermutedCoordinates coord(dim);
  enzo_float dtdx_i = dt/cell_width;

  auto accumulate = [dtdx_i,stale_depth,coord,
                     &flux_map,&dUcons_map](const std::string& key)
    {
      CSlice full_ax(nullptr, nullptr);
      // Since we don't have fluxes on the exterior faces along axis i, we can
      // not update dU in the first & last cell along the axis. Thus we cast
      // the calculation as:
      //   dU(k,j,i+1) -= dtdx_i * (flux(k, j, i+3/2) - flux(k, j, i+1/2))

      // define : dU_center(k,j,i) -> dU(k,j,i+1)
      EFlt3DArray dU, dU_center;
      dU = dUcons_map.get(key,stale_depth);
      dU_center = coord.get_subarray(dU, full_ax, full_ax, CSlice(1, -1));

      // define:  fl(k,j,i)        -> flux(k, j, i+1/2)
      //          fr(k,j,i)        -> flux(k, j, i+3/2)
      EFlt3DArray flux, fl, fr;
      flux = flux_map.get(key,stale_depth);
      fl = coord.get_subarray(flux, full_ax, full_ax, CSlice(0, -1));
      fr = coord.get_subarray(flux, full_ax, full_ax, CSlice(1, nullptr));

      for (int iz=0; iz<dU_center.shape(0); iz++) {
	for (int iy=0; iy<dU_center.shape(1); iy++) {
	  for (int ix=0; ix<dU_center.shape(2); ix++) {
	    dU_center(iz,iy,ix) -= dtdx_i * (fr(iz,iy,ix) - fl(iz,iy,ix));
	  }
	}
      }

    };

  for (const std::string& key : integration_keys_){ accumulate(key); }
  for (const std::string& key : passive_list){ accumulate(key); }
}

//----------------------------------------------------------------------

EFlt3DArray* EnzoIntegrationQuanUpdate::load_integration_quantities_
(EnzoEFltArrayMap &map, int stale_depth) const
{
  std::size_t nfields = integration_keys_.size();
  EFlt3DArray* arr = new EFlt3DArray[nfields];
  for (std::size_t i = 0; i<nfields; i++){
    arr[i] = map.get(integration_keys_[i], stale_depth);
  }
  return arr;
}

//----------------------------------------------------------------------

void EnzoIntegrationQuanUpdate::update_quantities
(EnzoEFltArrayMap &initial_integration_map, EnzoEFltArrayMap &dUcons_map,
 EnzoEFltArrayMap &out_integration_map,
 EnzoEquationOfState *eos, const int stale_depth,
 const str_vec_t &passive_list) const
{

  // Update passive scalars, it doesn't currently support renormalizing to 1
  update_passive_scalars_(initial_integration_map, dUcons_map,
                          out_integration_map, stale_depth, passive_list);

  // For now, not having density floor affect momentum or total energy density
  enzo_float density_floor = eos->get_density_floor();

  EFlt3DArray *cur_prim, *dU, *out_prim;
  cur_prim = load_integration_quantities_(initial_integration_map, stale_depth);
  dU       = load_integration_quantities_(dUcons_map, stale_depth);
  out_prim = load_integration_quantities_(out_integration_map,stale_depth);
  std::size_t nfields = integration_keys_.size();

  for (int iz = 1; iz < (cur_prim[density_index_].shape(0) - 1); iz++) {
    for (int iy = 1; iy < (cur_prim[density_index_].shape(1) - 1); iy++) {
      for (int ix = 1; ix < (cur_prim[density_index_].shape(2) - 1); ix++) {

	// get the initial density
	enzo_float old_rho = cur_prim[density_index_](iz,iy,ix);

	// now update the integration quantities that are conserved
	for (std::size_t i = 0; i < first_specific_index_; i++){
	  out_prim[i](iz,iy,ix) = cur_prim[i](iz,iy,ix) + dU[i](iz,iy,ix);
	}

	// possibly place a floor on the updated density.
	enzo_float new_rho = out_prim[density_index_](iz,iy,ix);
	new_rho = EnzoEquationOfState::apply_floor(new_rho, density_floor);
	out_prim[density_index_](iz,iy,ix) = new_rho;

	// update the specific integration primitives
	enzo_float inv_new_rho = 1./new_rho;
	for (std::size_t i = first_specific_index_; i < nfields; i++){
	  out_prim[i](iz,iy,ix) =
	    (cur_prim[i](iz,iy,ix) * old_rho + dU[i](iz,iy,ix)) * inv_new_rho;
	}

      }
    }
  }

  // apply floor to energy and sync the internal energy with total energy
  // (the latter only occurs if the dual energy formalism is in use)
  eos->apply_floor_to_energy_and_sync(out_integration_map, stale_depth + 1);

  delete[] cur_prim;  delete[] dU;  delete[] out_prim;
}

//----------------------------------------------------------------------

void EnzoIntegrationQuanUpdate::update_passive_scalars_
(EnzoEFltArrayMap &initial_integration_map, EnzoEFltArrayMap &dUcons_map,
 EnzoEFltArrayMap &out_integration_map, const int stale_depth,
 const str_vec_t &passive_list) const
{
  // cell-centered grid dimensions
  EFlt3DArray cur_rho = initial_integration_map.get("density", stale_depth);
  int mz, my, mx;
  mz = cur_rho.shape(0);    my = cur_rho.shape(1);    mx = cur_rho.shape(2);

  for (const std::string &key : passive_list){
    EFlt3DArray cur_conserved, out_conserved, dU;
    cur_conserved = initial_integration_map.get(key, stale_depth);
    out_conserved = out_integration_map.get(key, stale_depth);
    dU = dUcons_map.get(key, stale_depth);

    for (int iz=1; iz<mz-1; iz++) {
      for (int iy=1; iy<my-1; iy++) {
        for (int ix=1; ix<mx-1; ix++) {
          out_conserved(iz,iy,ix) = cur_conserved(iz,iy,ix) + dU(iz,iy,ix);
        }
      }
    }
  }

}
