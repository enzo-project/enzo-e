// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoIntegrableUpdate.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon June 24 2019
/// @brief    [\ref Enzo] Implementation of EnzoIntegrableUpdate.

#include "cello.hpp"
#include "enzo.hpp"

void append_grouping_pairs_(std::vector<std::string> integrable_groups,
			    FieldCat target_cat, std::size_t *density_index,
                            std::vector<std::string> &key_vec)
{
  int rank = cello::rank();
  for (std::string name : integrable_groups){
    bool vector_quantity, actively_advected;
    FieldCat category;
    bool success = EnzoCenteredFieldRegistry::quantity_properties
      (name, &vector_quantity, &category, &actively_advected);

    // Sanity Checks:
    ASSERT1("append_grouping_pairs_",
	    ("\"%s\" is not registered in EnzoCenteredFieldRegistry"),
	    name.c_str(), success);
    ASSERT1("append_grouping_pairs_",
	    ("\"%s\" should not be listed as an integrable quantity because "
	     "it is not actively advected."),
	    name.c_str(), actively_advected);

    if (category != target_cat){
      ASSERT1("append_grouping_pairs_",
	      ("Can't handle the integrable \"%s\" quantity because it has a "
	       "field category of FieldCat::other"),
	      name.c_str(), category != FieldCat::other);
      continue;
    } else if ((density_index != NULL) && (name == "density")){
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

EnzoIntegrableUpdate::EnzoIntegrableUpdate
(std::vector<std::string> integrable_groups, bool skip_B_update,
 std::vector<std::string> passive_groups) throw()
{
  passive_groups_ = passive_groups;

  if (skip_B_update){
    // remove bfield group from integrable_groups (if present)
    std::string target("bfield");
    integrable_groups.erase(std::remove(integrable_groups.begin(),
					integrable_groups.end(), target),
			    integrable_groups.end());
  }

  // assemble combined_integrable_groups_
  auto contains = [](const std::vector<std::string> &vec,
		     const std::string &value)
  { return std::find(vec.cbegin(),vec.cend(), value) != vec.cend(); };

  for (const std::string& elem : integrable_groups){
    if (!contains(combined_integrable_groups_, elem)) {
      combined_integrable_groups_.push_back(elem);
    }
  }
  for (const std::string& elem : passive_groups_){
    if (!contains(combined_integrable_groups_, elem)) {
      combined_integrable_groups_.push_back(elem);
    }
  }

  // sanity check:
  std::vector<std::string>::size_type num_unique_groups;
  num_unique_groups = this->combined_integrable_groups().size();
  ASSERT("EnzoIntegrableUpdate",
	 ("group names appear more than once in integrable_groups or "
	  "passive_groups"),
	 num_unique_groups == (integrable_groups.size() +
			       passive_groups_.size()) );

  // prepare integrable_keys_
  ASSERT("EnzoIntegrableUpdate",
	 ("\"density\" must be a registered integrable group."),
	 contains(integrable_groups, "density"));
  // First, add conserved quantities to integrable_keys_
  append_grouping_pairs_(integrable_groups, FieldCat::conserved,
			 &density_index_, integrable_keys_);
  first_specific_index_ = integrable_keys_.size();
  // Then, add specific quantities to integrable_keys_
  append_grouping_pairs_(integrable_groups, FieldCat::specific,
			 NULL, integrable_keys_);
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::clear_dUcons_map
(EnzoEFltArrayMap &dUcons_map, enzo_float value,
 const std::vector<std::vector<std::string>> &passive_lists) const noexcept
{
  auto clear_arr = [value,&dUcons_map](const std::string& key)
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

  for (const std::string& key : integrable_keys_){ clear_arr(key); }

  for (const std::vector<std::string>& cur_list : passive_lists){
    for (const std::string& key : cur_list){ clear_arr(key); }
  }
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::accumulate_flux_component
(int dim, double dt, double cell_width, EnzoEFltArrayMap &flux_map,
 EnzoEFltArrayMap &dUcons_map, int stale_depth,
 const std::vector<std::vector<std::string>> &passive_lists) const noexcept
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

  // first, handle actively integrated quantities:
  for (const std::string& key : integrable_keys_){ accumulate(key); }

  // second, handle passive scalars:
  for (const std::vector<std::string> key_list : passive_lists){
    for (const std::string& key : key_list){ accumulate(key); }
  }

}

//----------------------------------------------------------------------

EFlt3DArray* EnzoIntegrableUpdate::load_integrable_quantities_
(EnzoEFltArrayMap &map, int stale_depth) const
{
  std::size_t nfields = integrable_keys_.size();
  EFlt3DArray* arr = new EFlt3DArray[nfields];
  for (std::size_t i = 0; i<nfields; i++){
    arr[i] = map.get(integrable_keys_[i], stale_depth);
  }
  return arr;
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::update_quantities
(EnzoEFltArrayMap &initial_integrable_map, EnzoEFltArrayMap &dUcons_map,
 EnzoEFltArrayMap &out_integrable_map,
 EnzoEFltArrayMap &out_conserved_passive_scalar,
 EnzoEquationOfState *eos, int stale_depth,
 const std::vector<std::vector<std::string>> &passive_lists) const
{

  // Update passive scalars, it doesn't currently support renormalizing to 1
  update_passive_scalars_(initial_integrable_map, dUcons_map,
                          out_conserved_passive_scalar, stale_depth,
                          passive_lists);

  // For now, not having density floor affect momentum or total energy density
  enzo_float density_floor = eos->get_density_floor();

  EFlt3DArray *cur_prim, *dU, *out_prim;
  cur_prim = load_integrable_quantities_(initial_integrable_map, stale_depth);
  dU       = load_integrable_quantities_(dUcons_map, stale_depth);
  out_prim = load_integrable_quantities_(out_integrable_map,stale_depth);
  std::size_t nfields = integrable_keys_.size();

  for (int iz = 1; iz < (cur_prim[density_index_].shape(0) - 1); iz++) {
    for (int iy = 1; iy < (cur_prim[density_index_].shape(1) - 1); iy++) {
      for (int ix = 1; ix < (cur_prim[density_index_].shape(2) - 1); ix++) {

	// get the initial density
	enzo_float old_rho = cur_prim[density_index_](iz,iy,ix);

	// now update the integrable primitives that are conserved
	for (std::size_t i = 0; i < first_specific_index_; i++){
	  out_prim[i](iz,iy,ix) = cur_prim[i](iz,iy,ix) + dU[i](iz,iy,ix);
	}

	// possibly place a floor on the updated density.
	enzo_float new_rho = out_prim[density_index_](iz,iy,ix);
	new_rho = EnzoEquationOfState::apply_floor(new_rho, density_floor);
	out_prim[density_index_](iz,iy,ix) = new_rho;

	// update the specific integrable primitives
	enzo_float inv_new_rho = 1./new_rho;
	for (std::size_t i = first_specific_index_; i < nfields; i++){
	  out_prim[i](iz,iy,ix) =
	    (cur_prim[i](iz,iy,ix) * old_rho + dU[i](iz,iy,ix)) * inv_new_rho;
	}

      }
    }
  }

  delete[] cur_prim;  delete[] dU;  delete[] out_prim;
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::update_passive_scalars_
(EnzoEFltArrayMap &initial_integrable_map, EnzoEFltArrayMap &dUcons_map,
 EnzoEFltArrayMap &out_conserved_passive_scalar, int stale_depth,
 const std::vector<std::vector<std::string>> &passive_lists) const
{

  EFlt3DArray cur_rho = initial_integrable_map.get("density", stale_depth);

  // cell-centered grid dimensions
  int mz, my, mx;
  mz = cur_rho.shape(0);    my = cur_rho.shape(1);    mx = cur_rho.shape(2);

  std::vector<std::string> group_names = this->passive_groups_;
  for (std::size_t i = 0; i < passive_lists.size(); i++){

    const std::vector<std::string> &cur_list = passive_lists[i];

    if ((i > 0) && (cur_list.size() > 0)){
      ERROR("EnzoIntegrableUpdate::update_passive_scalars_",
            "This function does not currently support sets of passively "
            "advected scalars that must sum to 1.");
    }

    for (const std::string &key : cur_list){
      EFlt3DArray cur_specific, out_conserved, dU;
      cur_specific = initial_integrable_map.get(key, stale_depth);
      out_conserved = out_conserved_passive_scalar.get(key, stale_depth);
      dU = dUcons_map.get(key, stale_depth);

      for (int iz=1; iz<mz-1; iz++) {
	for (int iy=1; iy<my-1; iy++) {
	  for (int ix=1; ix<mx-1; ix++) {
	    out_conserved(iz,iy,ix)
	      = (cur_specific(iz,iy,ix) * cur_rho(iz,iy,ix) + dU(iz,iy,ix));
	  }
	}
      }

    }
  }
}
