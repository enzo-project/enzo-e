// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEOSIdeal.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of EnzoEOSIdeal

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoEOSIdeal::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change
  PUP::able::pup(p);
  p|gamma_;
  p|density_floor_;
  p|pressure_floor_;
}

//----------------------------------------------------------------------

// Helper function that performs a quick check to confirm that certain fields
// are  by both reconstructable_group and integrable_group
void confirm_same_fields_(Grouping &grouping_ref, Grouping &grouping_check,
			  std::string ref_name, std::string check_name,
			  std::vector<std::string> group_names,
			  std::string func_name)
{
  for (std::size_t i = 0; i<group_names.size(); i++){
    std::string group_name = group_names[i];
    int num_ref_fields = grouping_ref.size(group_name);
    int num_check_fields = grouping_check.size(group_name);

    ASSERT3(func_name.c_str(),
	    ("%s and %s groupings must have the same number of entries "
	     "for the %s group"),
	    ref_name.c_str(), check_name.c_str(), group_name.c_str(),
	    num_ref_fields == num_check_fields);

    for (int j=0;j<num_ref_fields;j++){
      std::string ref_field = grouping_ref.item(group_name,j);
      std::string check_field = grouping_check.item(group_name,j);

      ASSERT4(func_name.c_str(),
	      ("Field %d of the %s group is expected to have the same name in "
	       "the %s and %s groupings"),
	      j, group_name.c_str(), ref_name.c_str(), check_name.c_str(),
	      ref_field == check_field);
    }
  }
}

void check_recon_integ_overlap_(Grouping &reconstructable_group,
				Grouping &integrable_group,
				std::string func_name)
{
  // We assume that the following groups are represented by the same fields in
  // integrable and reconstructable
  std::vector<std::string> common_groups = {"density", "velocity", "bfield"};

  // we also expect overlap with the passive scalars
  EnzoCenteredFieldRegistry reg;
  std::vector<std::string> scalar_groups = reg.passive_scalar_group_names();
  common_groups.insert(common_groups.end(), scalar_groups.begin(),
		       scalar_groups.end());
    
  confirm_same_fields_(integrable_group, reconstructable_group,
		       "integrable", "reconstructable", common_groups,
		       func_name);
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::reconstructable_from_integrable
  (Block *block, Grouping &integrable_group, Grouping &reconstructable_group,
   Grouping &conserved_passive_group, int stale_depth)
{

  if (this->uses_dual_energy_formalism()){
    ERROR("EnzoEOSIdeal::reconstructable_from_integrable",
	  "dual energy formalism is not yet implemented.");
  }

  // Confirm that the expected fields (e.g. density, vx, vy, vz, bx, by, bz)
  // are the same in reconstructable_group and integrable_group
  check_recon_integ_overlap_(reconstructable_group, integrable_group,
			     "EnzoEOSIdeal::reconstructable_from_integrable");

  // Simply compute the pressure
  std::string pressure_name = reconstructable_group.item("pressure",0);
  pressure_from_integrable(block, integrable_group, pressure_name,
			   conserved_passive_group, stale_depth);
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::integrable_from_reconstructable
  (Block *block, Grouping &reconstructable_group, Grouping &integrable_group,
   int stale_depth, int reconstructed_axis)
{

  if (this->uses_dual_energy_formalism()){
    // If dual energy formalism is used then, then we need to also compute the
    // internal energy from the pressure

    // For now, we are just going to raise an error
    ERROR("EnzoEOSIdeal::integrable_from_reconstructable",
	  "dual energy formalism is not yet implemented.");
  }

  // Confirm that the expected fields (e.g. density, vx, vy, vz, bx, by, bz)
  // are the same in reconstructable_group and integrable_group 
  check_recon_integ_overlap_(reconstructable_group, integrable_group,
			     "EnzoEOSIdeal::integrable_from_reconstructable");

  EnzoFieldArrayFactory array_factory(block, stale_depth);

  // Define 2 temporary aliases to shorten code:
  Grouping recon_group = reconstructable_group;
  int rec_ax = reconstructed_axis;

  EFlt3DArray density, vx, vy, vz, pressure, bx, by, bz;
  density = retrieve_field_(array_factory, recon_group, "density", 0, rec_ax);
  vx = retrieve_field_(array_factory, recon_group, "velocity", 0, rec_ax);
  vy = retrieve_field_(array_factory, recon_group, "velocity", 1, rec_ax);
  vz = retrieve_field_(array_factory, recon_group, "velocity", 2, rec_ax);
  pressure = retrieve_field_(array_factory, recon_group, "pressure", 0,
			     rec_ax);
  bx = retrieve_field_(array_factory, recon_group, "bfield", 0, rec_ax);
  by = retrieve_field_(array_factory, recon_group, "bfield", 1, rec_ax);
  bz = retrieve_field_(array_factory, recon_group, "bfield", 2, rec_ax);

  EFlt3DArray etot = retrieve_field_(array_factory, integrable_group,
				     "total_energy", 0, rec_ax);
  
  ASSERT("EnzoEOSIdeal::integrable_from_reconstructable",
	 "Need to provide support for spatial variation in gamma",
	 !cello::field_descr()->is_field("gamma"));

  enzo_float inv_gm1 = 1./(get_gamma()-1.);

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	enzo_float v2 = (vx(iz,iy,ix) * vx(iz,iy,ix) +
			 vy(iz,iy,ix) * vy(iz,iy,ix) +
			 vz(iz,iy,ix) * vz(iz,iy,ix));
	enzo_float b2 = (bx(iz,iy,ix) * bx(iz,iy,ix) +
			 by(iz,iy,ix) * by(iz,iy,ix) +
			 bz(iz,iy,ix) * bz(iz,iy,ix));
	etot(iz,iy,ix) = ((pressure(iz,iy,ix) * inv_gm1
			   + 0.5 * b2) / density(iz,iy,ix)
			  + 0.5 * v2);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::pressure_from_integrable(Block *block,
					    Grouping &integrable_group,
					    std::string pressure_name,
					    Grouping &conserved_passive_group,
					    int stale_depth)
{

  // For now, we are not actually wrapping ComputePressure
  // To use EnzoComputePressure, we need to do some minor refactoring to allow
  // for optionally computing Pressure from fields specified in a Grouping.
  // This also requires making a modification to EnzoMethodGrackle's static
  // setup_grackle_fields method to also allow for specification of
  // relevant fields. Holding off on this for now
  //
  // As it stands, EnzoComputePressure (without Grackle) ALWAYS uses the
  // following fields to compute pressure 
  //   "density", "velocity_x", "velocity_y", "velocity_z", "total_energy",
  //   "bfield_x", "bfield_y", "bfield_z"
  
  if (enzo::config()->method_grackle_use_grackle){
    
    // I think internal energy may need to be precomputed (does that mean
    // we should always be using with dual energy formalism?) 
    ERROR("EnzoEOSIdeal::pressure_from_integrable",
	  "Not presently equipped to handle grackle");
  }

  if (this->uses_dual_energy_formalism()){
    ERROR("EnzoEOSIdeal::pressure_from_integrable",
  	  "dual energy formalism is not yet implemented.");
  }

  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EFlt3DArray density, vx, vy, vz, etot, bx, by, bz, pressure;

  // We are going to check that these are the specified fields
  density = array_factory.from_grouping(integrable_group, "density", 0);
  vx = array_factory.from_grouping(integrable_group, "velocity", 0);
  vy = array_factory.from_grouping(integrable_group, "velocity", 1);
  vz = array_factory.from_grouping(integrable_group, "velocity", 2);
  etot = array_factory.from_grouping(integrable_group, "total_energy", 0);
  bx = array_factory.from_grouping(integrable_group, "bfield", 0);
  by = array_factory.from_grouping(integrable_group, "bfield", 1);
  bz = array_factory.from_grouping(integrable_group, "bfield", 2);

  pressure = array_factory.from_name(pressure_name);

  enzo_float gm1 = get_gamma() - 1.;

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	enzo_float b2 = (bx(iz,iy,ix) * bx(iz,iy,ix) +
			 by(iz,iy,ix) * by(iz,iy,ix) +
			 bz(iz,iy,ix) * bz(iz,iy,ix));
	enzo_float v2 = (vx(iz,iy,ix) * vx(iz,iy,ix) +
			 vy(iz,iy,ix) * vy(iz,iy,ix) +
			 vz(iz,iy,ix) * vz(iz,iy,ix));
	pressure(iz,iy,ix) = 
	  gm1 * ((etot(iz,iy,ix) - 0.5 * v2) * density(iz,iy,ix) - 0.5 * b2);
      }
    }
  }
}


void EnzoEOSIdeal::pressure_from_reconstructable
(Block *block, Grouping &reconstructable_group, std::string pressure_name,
 int stale_depth, int reconstructed_axis)
{
  // This is necessary since other equations of state may not include pressure
  // as a reconstructable quantity.

  // We will check if pressure_name the field containing pressure in
  // reconstructable_group are the same, if so then do nothing. Otherwise,
  // simply copy the values over.

  if (reconstructable_group.item("pressure",0) == pressure_name){
    // The fields have the same name, we don't have to do anything
    return;
  } else {
    EnzoFieldArrayFactory array_factory(block, stale_depth);
    EFlt3DArray old_p = retrieve_field_(array_factory, reconstructable_group,
					"pressure", 0, reconstructed_axis);
    
    Grouping temp_group;
    temp_group.add(pressure_name, "pressure");
    EFlt3DArray new_p = retrieve_field_(array_factory, temp_group, "pressure",
					0, reconstructed_axis);

    for (int iz=0; iz< old_p.shape(0); iz++) {
      for (int iy=0; iy< old_p.shape(1); iy++) {
	for (int ix=0; ix< old_p.shape(2); ix++) {
	  new_p(iz,iy,ix) = old_p(iz,iy,ix);
	}
      }
    }
    // new_p.subarray() = old_p;
  }
}


//----------------------------------------------------------------------

// Applies the pressure_floor to total_energy
void EnzoEOSIdeal::apply_floor_to_total_energy(Block *block,
					       Grouping &integrable_group,
					       int stale_depth)
{

  if (enzo::config()->method_grackle_use_grackle){
    ERROR("EnzoEOSIdeal::apply_floor_to_total_energy",
	  "Not presently equipped to handle grackle");
    // since we are not currently allowing Grackle (due to possibly
    // variable gamma)
  }

  EnzoFieldArrayFactory array_factory(block,stale_depth);
  EFlt3DArray density, vx, vy, vz, etot, bx, by, bz;
  density = array_factory.from_grouping(integrable_group, "density", 0);
  vx = array_factory.from_grouping(integrable_group, "velocity", 0);
  vy = array_factory.from_grouping(integrable_group, "velocity", 1);
  vz = array_factory.from_grouping(integrable_group, "velocity", 2);
  etot = array_factory.from_grouping(integrable_group, "total_energy", 0);
  bx = array_factory.from_grouping(integrable_group, "bfield", 0);
  by = array_factory.from_grouping(integrable_group, "bfield", 1);
  bz = array_factory.from_grouping(integrable_group, "bfield", 2);

  enzo_float pressure = get_pressure_floor();
  enzo_float inv_gm1 = 1./(get_gamma()-1.);

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	enzo_float kinetic, magnetic;
	kinetic = 0.5*(vx(iz,iy,ix) * vx(iz,iy,ix) +
		       vy(iz,iy,ix) * vy(iz,iy,ix) +
		       vz(iz,iy,ix) * vz(iz,iy,ix));
	magnetic = 0.5*(bx(iz,iy,ix) * bx(iz,iy,ix) +
			by(iz,iy,ix) * by(iz,iy,ix) +
			bz(iz,iy,ix) * bz(iz,iy,ix))/density(iz,iy,ix);
	etot(iz,iy,ix) = std::max((pressure * inv_gm1/density(iz,iy,ix)
				   + kinetic + magnetic),
				  etot(iz,iy,ix));
      }
    }
  }
}

//----------------------------------------------------------------------

EFlt3DArray EnzoEOSIdeal::retrieve_field_(EnzoFieldArrayFactory &array_factory,
					  Grouping &group,
					  std::string group_name, int index,
					  int reconstructed_axis)
{
  if (reconstructed_axis == -1){
    return array_factory.from_grouping(group, group_name, index);
  } else {
    return array_factory.reconstructed_field(group, group_name, index,
					     reconstructed_axis);
  }
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::copy_passively_advected_fields_
(EnzoFieldArrayFactory &array_factory, Grouping &origin_group,
 Grouping &destination_group, int reconstruction_axis)
{
  EnzoCenteredFieldRegistry registry;
  std::vector<std::string> group_names = registry.passive_scalar_group_names();
  for (unsigned int i=0;i<group_names.size();i++){
    std::string group_name = group_names[i];
    int num_fields = origin_group.size(group_name);
    for (int j=0;j<num_fields;j++){

      EFlt3DArray src = retrieve_field_(array_factory, origin_group,
					group_name, j, reconstruction_axis);
      EFlt3DArray dest = retrieve_field_(array_factory, destination_group,
					 group_name, j, reconstruction_axis);
      dest.subarray() = src;
    }
  }
}

