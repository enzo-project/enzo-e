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

// Would eventually like this to just wrap the compute_pressure object
void EnzoEOSIdeal::compute_pressure(Block *block, Grouping &cons_group,
				    Grouping &prim_group, int stale_depth)
{
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EFlt3DArray density, p_x, p_y, p_z, etot, b_x, b_y, b_z;
  density = array_factory.from_grouping(cons_group, "density", 0);
  p_x = array_factory.from_grouping(cons_group, "momentum", 0);
  p_y = array_factory.from_grouping(cons_group, "momentum", 1);
  p_z = array_factory.from_grouping(cons_group, "momentum", 2);
  etot = array_factory.from_grouping(cons_group, "total_energy", 0);
  b_x = array_factory.from_grouping(cons_group, "bfield", 0);
  b_y = array_factory.from_grouping(cons_group, "bfield", 1);
  b_z = array_factory.from_grouping(cons_group, "bfield", 2);

  EFlt3DArray pressure;
  pressure = array_factory.from_grouping(prim_group, "pressure", 0);

  enzo_float gm1 = get_gamma() - 1.;
  // No good way to tell if fields are face-centered
  // Since, fields are the same size, this just means that some unnecessary
  // values are calculated
  // Iteration limits compatible with both 2D and 3D grids
  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	enzo_float magnetic, kinetic;
	magnetic = 0.5*(b_x(iz,iy,ix) * b_x(iz,iy,ix) +
			b_y(iz,iy,ix) * b_y(iz,iy,ix) +
			b_z(iz,iy,ix) * b_z(iz,iy,ix));
	kinetic = 0.5*(p_x(iz,iy,ix) * p_x(iz,iy,ix) +
		       p_y(iz,iy,ix) * p_y(iz,iy,ix) +
		       p_z(iz,iy,ix) * p_z(iz,iy,ix))/density(iz,iy,ix);
	pressure(iz,iy,ix) = gm1 * (etot(iz,iy,ix) - magnetic - kinetic);

      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::primitive_from_conservative(Block *block,
					       Grouping &cons_group,
					       Grouping &prim_group,
					       int stale_depth)
{
  compute_pressure(block, cons_group, prim_group, stale_depth);
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EFlt3DArray cons_density, px, py, pz, cons_bx, cons_by, cons_bz;
  cons_density = array_factory.from_grouping(cons_group, "density", 0);
  px = array_factory.from_grouping(cons_group, "momentum", 0);
  py = array_factory.from_grouping(cons_group, "momentum", 1);
  pz = array_factory.from_grouping(cons_group, "momentum", 2);
  cons_bx = array_factory.from_grouping(cons_group, "bfield", 0);
  cons_by = array_factory.from_grouping(cons_group, "bfield", 1);
  cons_bz = array_factory.from_grouping(cons_group, "bfield", 2);

  EFlt3DArray prim_density, vx, vy, vz, prim_bx, prim_by, prim_bz;
  prim_density = array_factory.from_grouping(prim_group, "density", 0);
  vx = array_factory.from_grouping(prim_group, "velocity", 0);
  vy = array_factory.from_grouping(prim_group, "velocity", 1);
  vz = array_factory.from_grouping(prim_group, "velocity", 2);
  prim_bx = array_factory.from_grouping(prim_group, "bfield", 0);
  prim_by = array_factory.from_grouping(prim_group, "bfield", 1);
  prim_bz = array_factory.from_grouping(prim_group, "bfield", 2);

  for (int iz=0; iz<cons_density.shape(0); iz++) {
    for (int iy=0; iy<cons_density.shape(1); iy++) {
      for (int ix=0; ix<cons_density.shape(2); ix++) {

	enzo_float density = cons_density(iz,iy,ix);
	prim_density(iz,iy,ix) = density;
	vx(iz,iy,ix) = px(iz,iy,ix)/density;
	vy(iz,iy,ix) = py(iz,iy,ix)/density;
	vz(iz,iy,ix) = pz(iz,iy,ix)/density;
	prim_bx(iz,iy,ix) = cons_bx(iz,iy,ix);
	prim_by(iz,iy,ix) = cons_by(iz,iy,ix);
	prim_bz(iz,iy,ix) = cons_bz(iz,iy,ix);
      }
    }
  }
  // Copy conserved passive fields to primitive
  copy_passively_advected_fields_(array_factory, cons_group, prim_group);
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::conservative_from_primitive(Block *block,
					       Grouping &prim_group,
					       Grouping &cons_group,
					       int stale_depth,
					       int reconstructed_axis)
{
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EFlt3DArray prim_density, vx, vy, vz, pressure;
  EFlt3DArray prim_bx, prim_by, prim_bz;

  int rec_ax = reconstructed_axis;
  prim_density = retrieve_field_(array_factory, prim_group, "density", 0,
				 rec_ax);
  vx = retrieve_field_(array_factory, prim_group, "velocity", 0, rec_ax);
  vy = retrieve_field_(array_factory, prim_group, "velocity", 1, rec_ax);
  vz = retrieve_field_(array_factory, prim_group, "velocity", 2, rec_ax);
  pressure = retrieve_field_(array_factory, prim_group, "pressure", 0, rec_ax);
  prim_bx = retrieve_field_(array_factory, prim_group, "bfield", 0, rec_ax);
  prim_by = retrieve_field_(array_factory, prim_group, "bfield", 1, rec_ax);
  prim_bz = retrieve_field_(array_factory, prim_group, "bfield", 2, rec_ax);

  EFlt3DArray cons_density, px, py, pz, etot;
  EFlt3DArray cons_bx, cons_by, cons_bz;
  cons_density = retrieve_field_(array_factory, cons_group, "density", 0,
				 rec_ax);
  px = retrieve_field_(array_factory, cons_group, "momentum", 0, rec_ax);
  py = retrieve_field_(array_factory, cons_group, "momentum", 1, rec_ax);
  pz = retrieve_field_(array_factory, cons_group, "momentum", 2, rec_ax);
  etot = retrieve_field_(array_factory, cons_group, "total_energy", 0, rec_ax);
  cons_bx = retrieve_field_(array_factory, cons_group, "bfield", 0, rec_ax);
  cons_by = retrieve_field_(array_factory, cons_group, "bfield", 1, rec_ax);
  cons_bz = retrieve_field_(array_factory, cons_group, "bfield", 2, rec_ax);

  enzo_float inv_gm1 = 1./(get_gamma()-1.);

  // It's okay that this function is used for computing face-centered temporary
  // interior fields
  for (int iz=0; iz<cons_density.shape(0); iz++) {
    for (int iy=0; iy<cons_density.shape(1); iy++) {
      for (int ix=0; ix<cons_density.shape(2); ix++) {

	enzo_float density = prim_density(iz,iy,ix);

	cons_density(iz,iy,ix) = density;
	px(iz,iy,ix) = vx(iz,iy,ix)*density;
	py(iz,iy,ix) = vy(iz,iy,ix)*density;
	pz(iz,iy,ix) = vz(iz,iy,ix)*density;

	enzo_float kinetic, magnetic;
	kinetic = 0.5*(vx(iz,iy,ix) * vx(iz,iy,ix) +
		       vy(iz,iy,ix) * vy(iz,iy,ix) +
		       vz(iz,iy,ix) * vz(iz,iy,ix))*density;
	magnetic = 0.5*(prim_bx(iz,iy,ix) * prim_bx(iz,iy,ix) +
			prim_by(iz,iy,ix) * prim_by(iz,iy,ix) +
			prim_bz(iz,iy,ix) * prim_bz(iz,iy,ix));
	etot(iz,iy,ix) = pressure(iz,iy,ix) * inv_gm1 + kinetic + magnetic;

	cons_bx(iz,iy,ix) = prim_bx(iz,iy,ix);
	cons_by(iz,iy,ix) = prim_by(iz,iy,ix);
	cons_bz(iz,iy,ix) = prim_bz(iz,iy,ix);
      }
    }
  }
  // Copy primitive passive fields to conserved
  copy_passively_advected_fields_(array_factory, prim_group, cons_group,
				  reconstructed_axis);
}

//----------------------------------------------------------------------

// Applies the pressure_floor to total_energy
void EnzoEOSIdeal::apply_floor_to_energy(Block *block, Grouping &cons_group,
					 int stale_depth)
{
  EnzoFieldArrayFactory array_factory(block,stale_depth);
  EFlt3DArray density, px, py, pz, etot, bx, by, bz;
  density = array_factory.from_grouping(cons_group, "density", 0);
  px = array_factory.from_grouping(cons_group, "momentum", 0);
  py = array_factory.from_grouping(cons_group, "momentum", 1);
  pz = array_factory.from_grouping(cons_group, "momentum", 2);
  etot = array_factory.from_grouping(cons_group, "total_energy", 0);
  bx = array_factory.from_grouping(cons_group, "bfield", 0);
  by = array_factory.from_grouping(cons_group, "bfield", 1);
  bz = array_factory.from_grouping(cons_group, "bfield", 2);

  enzo_float pressure = get_pressure_floor();
  enzo_float inv_gm1 = 1./(get_gamma()-1.);
  
  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	enzo_float kinetic, magnetic;
	kinetic = 0.5*(px(iz,iy,ix) * px(iz,iy,ix) +
		       py(iz,iy,ix) * py(iz,iy,ix) +
		       pz(iz,iy,ix) * pz(iz,iy,ix))/density(iz,iy,ix);
	magnetic = 0.5*(bx(iz,iy,ix) * bx(iz,iy,ix) +
			by(iz,iy,ix) * by(iz,iy,ix) +
			bz(iz,iy,ix) * bz(iz,iy,ix));
	etot(iz,iy,ix) = std::max(pressure * inv_gm1 + kinetic + magnetic,
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

