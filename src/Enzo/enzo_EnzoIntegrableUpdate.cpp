// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoIntegrableUpdate.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon June 24 2019
/// @brief    [\ref Enzo] Implementation of EnzoIntegrableUpdate.

#include "cello.hpp"
#include "enzo.hpp"

void append_grouping_pairs_(std::vector<std::string> integrable_groups,
			    FieldCat target_cat, std::size_t *density_index,
			    EnzoCenteredFieldRegistry &registry,
			    std::vector<std::pair<std::string,int>> &pair_vec)
{
  int rank = cello::rank();
  for (std::string name : integrable_groups){
    bool vector_quantity, actively_advected;
    FieldCat category;
    registry.quantity_properties(name, &vector_quantity, &category,
				 &actively_advected);
    // Sanity Check
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
      *density_index = pair_vec.size();
    }

    int nfields = vector_quantity ? 3 : 1;
    for (int i = 0; i < nfields; i++){
      std::pair<std::string,int> pair(name,i);
      pair_vec.push_back(pair);
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

  // prepare integrable_grouping_items_
  ASSERT("EnzoIntegrableUpdate",
	 ("\"density\" must be a registered integrable group."),
	 contains(integrable_groups, "density"));
  // First, add conserved quantities to integrable_grouping_items_
  EnzoCenteredFieldRegistry reg;
  append_grouping_pairs_(integrable_groups, FieldCat::conserved,
			 &density_index_, reg, integrable_grouping_items_);
  first_specific_index_ = integrable_grouping_items_.size();
  append_grouping_pairs_(integrable_groups, FieldCat::specific,
			 NULL, reg, integrable_grouping_items_);
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::clear_dUcons_group(Block *block,
					      Grouping &dUcons_group,
					      enzo_float value) const
{
  const std::vector<std::string> group_names =combined_integrable_groups();
  EnzoFieldArrayFactory array_factory(block, 0);

  for (const std::string& name : group_names){
    int num_fields = dUcons_group.size(name);

    for (int i=0; i<num_fields; i++){
      EFlt3DArray array = array_factory.from_grouping(dUcons_group, name, i);

      for (int iz=0; iz<array.shape(0); iz++) {
	for (int iy=0; iy<array.shape(1); iy++) {
	  for (int ix=0; ix<array.shape(2); ix++) {
	    array(iz,iy,ix) = value;
	  }
	}
      }

    }
  }
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::accumulate_flux_component(Block *block,
						     int dim, double dt,
						     Grouping &flux_group,
						     Grouping &dUcons_group,
						     int stale_depth) const
{
  const std::vector<std::string> group_names =combined_integrable_groups();
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EnzoPermutedCoordinates coord(dim);

  EnzoBlock * enzo_block = enzo::block(block);
  enzo_float dtdx_i = dt/enzo_block->CellWidth[coord.i_axis()];

  CSlice full_ax(nullptr, nullptr);

  for (std::string name : group_names){

    int num_fields = dUcons_group.size(name);
    ASSERT1("EnzoIntegrableUpdate::accumulate_flux_divergence",
	    ("The number of fields in the \"%s\" group is differnt between "
	     "flux_group and dUcons_group."), name.c_str(),
	    num_fields == flux_group.size(name));

    for (int field_ind = 0; field_ind < num_fields; field_ind++){
      // Since we don't have fluxes on the exterior faces along axis i, we can
      // not update dU in the first & last cell along the axis. Thus we cast
      // the calculation as:
      //   dU(k,j,i+1) -= dtdx_i * (flux(k, j, i+3/2) - flux(k, j, i+1/2))

      // define : dU_center(k,j,i) -> dU(k,j,i+1)
      EFlt3DArray dU, dU_center;
      dU = array_factory.from_grouping(dUcons_group, name, field_ind);
      dU_center = coord.get_subarray(dU, full_ax, full_ax, CSlice(1, -1));

      // define:  fl(k,j,i)        -> flux(k, j, i+1/2)
      //          fr(k,j,i)        -> flux(k, j, i+3/2)
      EFlt3DArray flux, fl, fr;
      flux = array_factory.from_grouping(flux_group, name, field_ind);
      fl = coord.get_subarray(flux, full_ax, full_ax, CSlice(0, -1));
      fr = coord.get_subarray(flux, full_ax, full_ax, CSlice(1, nullptr));

      for (int iz=0; iz<dU_center.shape(0); iz++) {
	for (int iy=0; iy<dU_center.shape(1); iy++) {
	  for (int ix=0; ix<dU_center.shape(2); ix++) {
	    dU_center(iz,iy,ix) -= dtdx_i * (fr(iz,iy,ix) - fl(iz,iy,ix));
	  }
	}
      }

    }
  }
}

//----------------------------------------------------------------------

EFlt3DArray* EnzoIntegrableUpdate::load_integrable_quantites_
(Block *block, Grouping & grouping, int stale_depth) const
{
  std::size_t nfields = integrable_grouping_items_.size();
  EFlt3DArray* arr = new EFlt3DArray[nfields];
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  for (std::size_t i = 0; i<nfields; i++){
    std::pair<std::string, int> pair = integrable_grouping_items_[i];
    arr[i] = array_factory.from_grouping(grouping, pair.first, pair.second);
  }
  return arr;
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::update_quantities
(Block *block, Grouping &initial_integrable_group, Grouping &dUcons_group,
 Grouping &out_integrable_group, Grouping &out_conserved_passive_scalar,
 EnzoEquationOfState *eos, int stale_depth) const
{

  // Update passive scalars, it doesn't currently support renormalizing to 1
  update_passive_scalars_(block, initial_integrable_group, dUcons_group,
			  out_conserved_passive_scalar, stale_depth);

  // For now, not having density floor affect momentum or total energy density
  enzo_float density_floor = eos->get_density_floor();

  EnzoCenteredFieldRegistry registry;
  EFlt3DArray *cur_prim, *dU, *out_prim;
  cur_prim = load_integrable_quantites_(block, initial_integrable_group,
					stale_depth);
  dU       = load_integrable_quantites_(block, dUcons_group, stale_depth);
  out_prim = load_integrable_quantites_(block, out_integrable_group,
					stale_depth);
  std::size_t nfields = integrable_grouping_items_.size();

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

  // apply floor to energy and sync the internal energy with total energy
  // (the latter only occurs if the dual energy formalism is in use)
  eos->apply_floor_to_energy_and_sync(block, out_integrable_group,
				      stale_depth+1);

  delete[] cur_prim;  delete[] dU;  delete[] out_prim;
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::update_passive_scalars_
  (Block *block, Grouping &initial_integrable_group, Grouping &dUcons_group,
   Grouping &out_conserved_passive_scalar, int stale_depth) const
{

  // Does not currently handle passively advected scalars that must sum to 1
  // The easiest way to do that would be to divide each group of fields that
  // must sum to 1 into different groupings

  EnzoFieldArrayFactory array_factory(block, stale_depth);

  EFlt3DArray cur_rho
    = array_factory.from_grouping(initial_integrable_group, "density", 0);

  // cell-centered grid dimensions
  int mz, my, mx;
  mz = cur_rho.shape(0);    my = cur_rho.shape(1);    mx = cur_rho.shape(2);

  std::vector<std::string> group_names = this->passive_groups_;
  for (std::string group_name : group_names){
    int num_fields = initial_integrable_group.size(group_name);

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      EFlt3DArray cur_specific, out_conserved, dU;
      cur_specific = array_factory.from_grouping(initial_integrable_group,
						 group_name, field_ind);
      out_conserved = array_factory.from_grouping(out_conserved_passive_scalar,
						  group_name, field_ind);
      dU = array_factory.from_grouping(dUcons_group, group_name, field_ind);

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
