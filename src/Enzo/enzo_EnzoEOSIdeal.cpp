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
   int stale_depth)
{

  if (this->uses_dual_energy_formalism()){
    // If dual energy formalism is used then nothing should happen here
    // Just to be safe though, we are raising this error for now
    ERROR("EnzoEOSIdeal::reconstructable_from_integrable",
	  "dual energy formalism is not yet implemented.");

    // In this case we should confirm that the same internal energy field is
    // tracked in both cases
  }

  // Confirm that the relevant fields overlap as expected
  check_recon_integ_overlap_(reconstructable_group, integrable_group,
			     "EnzoEOSIdeal::reconstructable_from_integrable");

  EnzoFieldArrayFactory array_factory(block, stale_depth);

  EFlt3DArray density, vx, vy, vz, etot, bx, by, bz;
  density = array_factory.from_grouping(integrable_group, "density", 0);
  vx = array_factory.from_grouping(integrable_group, "velocity", 0);
  vy = array_factory.from_grouping(integrable_group, "velocity", 1);
  vz = array_factory.from_grouping(integrable_group, "velocity", 2);
  etot = array_factory.from_grouping(integrable_group, "total_energy", 0);
  bx = array_factory.from_grouping(integrable_group, "bfield", 0);
  by = array_factory.from_grouping(integrable_group, "bfield", 1);
  bz = array_factory.from_grouping(integrable_group, "bfield", 2);

  EFlt3DArray eint = array_factory.from_grouping(reconstructable_group,
						 "internal_energy", 0);

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	// compute specific kinetic and magnetic energies
	enzo_float kinetic, magnetic;
	kinetic = 0.5*(vx(iz,iy,ix) * vx(iz,iy,ix) +
		       vy(iz,iy,ix) * vy(iz,iy,ix) +
		       vz(iz,iy,ix) * vz(iz,iy,ix));
	magnetic = 0.5*(bx(iz,iy,ix) * bx(iz,iy,ix) +
			by(iz,iy,ix) * by(iz,iy,ix) +
			bz(iz,iy,ix) * bz(iz,iy,ix)) / density(iz,iy,ix);
	eint(iz,iy,ix) = etot(iz,iy,ix) - kinetic - magnetic;
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::integrable_from_reconstructable
  (Block *block, Grouping &reconstructable_group, Grouping &integrable_group,
   int stale_depth, int reconstructed_axis)
{

  if (this->uses_dual_energy_formalism()){
    // If dual energy formalism is used then nothing different should happen
    // here. We still need to compute total energy from specific internal
    // energy
    // However, we should confirm that the internal energy field matches for
    // both fields

    // For now, we are just going to raise an error
    ERROR("EnzoEOSIdeal::integrable_from_reconstructable",
	  "dual energy formalism is not yet implemented.");
  }


  // Confirm that the relevant fields overlap as expected
  check_recon_integ_overlap_(reconstructable_group, integrable_group,
			     "EnzoEOSIdeal::integrable_from_reconstructable");

  EnzoFieldArrayFactory array_factory(block, stale_depth);

  // Define 2 temporary aliases to shorten code:
  Grouping recon_group = reconstructable_group;
  int rec_ax = reconstructed_axis;

  EFlt3DArray density, vx, vy, vz, eint, bx, by, bz;

  density = retrieve_field_(array_factory, recon_group, "density", 0, rec_ax);
  vx = retrieve_field_(array_factory, recon_group, "velocity", 0, rec_ax);
  vy = retrieve_field_(array_factory, recon_group, "velocity", 1, rec_ax);
  vz = retrieve_field_(array_factory, recon_group, "velocity", 2, rec_ax);
  eint = retrieve_field_(array_factory, recon_group, "internal_energy", 0,
			 rec_ax);
  bx = retrieve_field_(array_factory, recon_group, "bfield", 0, rec_ax);
  by = retrieve_field_(array_factory, recon_group, "bfield", 1, rec_ax);
  bz = retrieve_field_(array_factory, recon_group, "bfield", 2, rec_ax);

  EFlt3DArray etot = retrieve_field_(array_factory, integrable_group,
				     "total_energy", 0, rec_ax);

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	// compute specific kinetic and magnetic energies
	enzo_float kinetic, magnetic;
	kinetic = 0.5*(vx(iz,iy,ix) * vx(iz,iy,ix) +
		       vy(iz,iy,ix) * vy(iz,iy,ix) +
		       vz(iz,iy,ix) * vz(iz,iy,ix));
	magnetic = 0.5*(bx(iz,iy,ix) * bx(iz,iy,ix) +
			by(iz,iy,ix) * by(iz,iy,ix) +
			bz(iz,iy,ix) * bz(iz,iy,ix)) / density(iz,iy,ix);
	etot(iz,iy,ix) = eint(iz,iy,ix) + kinetic + magnetic;
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::pressure_from_integrable(Block *block,
					    Grouping &integrable_group,
					    std::string pressure_name,
					    Grouping &passive_scalars_group,
					    bool specific_passive_scalars,
					    int stale_depth)
{

  if (specific_passive_scalars){
    // I am currently unaware of a scenario where we might need to compute
    // thermal pressure from specific_passive_scalars, but I wanted to add it to
    // the function signature to be explicit
    ERROR("EnzoEOSIdeal::pressure_from_integrable",
	  "Not presently equipped to compute pressure with specific scalars");
  }

  if (this->uses_dual_energy_formalism()){
    ERROR("EnzoEOSIdeal::pressure_from_integrable",
	  "dual energy formalism is not yet implemented.");
  }
  
  if (enzo::config()->method_grackle_use_grackle){
    // I think internal energy may need to be precomputed (does that mean
    // we should always be using with dual energy formalism?) 
    ERROR("EnzoEOSIdeal::pressure_from_integrable",
	  "Not presently equipped to handle grackle");

    // since we are not currently allowing Grackle,  passive scalars are
    // unimportant
  }

  // For now, this method just forwards to EnzoComputePressure without
  // allowing Grackle. In this case EnzoComputePressure only uses the fields
  //   "density", "velocity_x", "velocity_y", "velocity_z", "total_energy",
  //   "bfield_x", "bfield_y", "bfield_z"

  // We are going to check that these are the specified fields

  ASSERT("EnzoEOSIdeal::pressure_from_integrable",
	 "Currently, the supplied grouping must include \"density\"",
	 integrable_group.is_in( "density", "density"));
  ASSERT("EnzoEOSIdeal::pressure_from_integrable",
	 "Currently, the supplied grouping must include \"total_energy\"",
	 integrable_group.is_in( "total_energy", "total_energy"));

  std::vector<std::string> vel{"velocity_x", "velocity_y", "velocity_z"};
  std::vector<std::string> b{"bfield_x", "bfield_y", "bfield_z"};

  for (std::size_t i=0; i<3; i++){
    std::string velocity_name = vel[i];
    std::string bfield_name   = b[i];
    ASSERT1("EnzoEOSIdeal::pressure_from_integrable",
	    "Currently, the supplied grouping must include \"%s\"",
	    vel[i].c_str(), integrable_group.is_in( vel[i], "velocity"));
    ASSERT1("EnzoEOSIdeal::pressure_from_integrable",
	    "Currently, the supplied grouping must include \"%s\"",
	    b[i].c_str(), integrable_group.is_in( b[i], "bfield"));
  }

  // assumes that we are not using comoving coordinates
  EnzoComputePressure compute_pressure (this->get_gamma(),false);
  Field field = block->data()->field();
  enzo_float* pressure = (enzo_float*) field.values(pressure_name);
  compute_pressure.compute(block, pressure);
}


void EnzoEOSIdeal::pressure_from_reconstructable
(Block *block, Grouping &reconstructable_group, std::string pressure_name,
 int stale_depth, int reconstructed_axis)
{
  // Since the reconstructable primitives never include total_energy (just
  // internal energy), there is no difference when the dual energy formalism is
  // adopted

  // A modification may have to be made from having a variable specific heat
  // ratio

  EnzoFieldArrayFactory array_factory(block, stale_depth);
  // Load the reconstructed values
  EFlt3DArray density, eint;
  density = retrieve_field_(array_factory, reconstructable_group, "density",
			    0, reconstructed_axis);
  eint = retrieve_field_(array_factory, reconstructable_group,
			 "internal_energy",  0, reconstructed_axis);

  // load the pressure - the way we do this is kind of silly (defining a
  // temporary grouping since reconstructable_group is probably reconstructed), 
  // but not currently worth modifying the interface of EnzoFieldArrayFactory
  Grouping temp_group;
  temp_group.add(pressure_name, "pressure");
  EFlt3DArray pressure = retrieve_field_(array_factory, temp_group, "pressure",
					 0, reconstructed_axis);

  enzo_float gm1 = (get_gamma()-1.);

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	pressure(iz,iy,ix) = gm1 * eint(iz,iy,ix) * density(iz, iy, ix);

      }
    }
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

// Applies the pressure_floor to internal_energy
void EnzoEOSIdeal::apply_floor_to_internal_energy
(Block *block, Grouping &reconstructable_group, int stale_depth)
{

  if (enzo::config()->method_grackle_use_grackle){
    ERROR("EnzoEOSIdeal::apply_floor_to_internal_energy",
	  "Not presently equipped to handle grackle");
    // since we are not currently allowing Grackle (due to possibly
    // variable gamma)
  }
  
  EnzoFieldArrayFactory array_factory(block,stale_depth);
  EFlt3DArray density = array_factory.from_grouping(reconstructable_group,
						    "density", 0);
  EFlt3DArray eint = array_factory.from_grouping(reconstructable_group,
						 "internal_energy", 0);
  enzo_float pressure = get_pressure_floor();
  enzo_float inv_gm1 = 1./(get_gamma()-1.);

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	enzo_float floor = pressure*inv_gm1/density(iz,iy,ix);

	eint(iz,iy,ix) = EnzoEquationOfState::apply_floor(eint(iz,iy,ix),
							  floor);
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

