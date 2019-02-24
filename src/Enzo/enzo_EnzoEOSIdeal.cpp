
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
void EnzoEOSIdeal::compute_pressure(Block *block,
				    Grouping &cons_group,
				    Grouping &prim_group)
{
  EnzoFieldArrayFactory array_factory(block);
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
					       Grouping &prim_group)
{
  compute_pressure(block, cons_group, prim_group);
  EnzoFieldArrayFactory array_factory(block);
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
}

//----------------------------------------------------------------------


void EnzoEOSIdeal::conservative_from_primitive(Block *block,
					       Grouping &prim_group,
					       Grouping &cons_group)
{
  EnzoFieldArrayFactory array_factory(block);
  EFlt3DArray prim_density, vx, vy, vz, pressure;
  EFlt3DArray prim_bx, prim_by, prim_bz;
  prim_density = array_factory.from_grouping(prim_group, "density", 0);
  vx = array_factory.from_grouping(prim_group, "velocity", 0);
  vy = array_factory.from_grouping(prim_group, "velocity", 1);
  vz = array_factory.from_grouping(prim_group, "velocity", 2);
  pressure = array_factory.from_grouping(prim_group, "pressure", 0);
  prim_bx = array_factory.from_grouping(prim_group, "bfield", 0);
  prim_by = array_factory.from_grouping(prim_group, "bfield", 1);
  prim_bz = array_factory.from_grouping(prim_group, "bfield", 2);

  EFlt3DArray cons_density, px, py, pz, etot;
  EFlt3DArray cons_bx, cons_by, cons_bz;
  cons_density = array_factory.from_grouping(cons_group, "density", 0);
  px = array_factory.from_grouping(cons_group, "momentum", 0);
  py = array_factory.from_grouping(cons_group, "momentum", 1);
  pz = array_factory.from_grouping(cons_group, "momentum", 2);
  etot = array_factory.from_grouping(cons_group, "total_energy", 0);
  cons_bx = array_factory.from_grouping(cons_group, "bfield", 0);
  cons_by = array_factory.from_grouping(cons_group, "bfield", 1);
  cons_bz = array_factory.from_grouping(cons_group, "bfield", 2);

  enzo_float inv_gm1 = 1./(get_gamma()-1.);

  // It should be okay that this function is being primarily used for
  // face-centered temporary interior fields
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
}

//----------------------------------------------------------------------

enzo_float EnzoEOSIdeal::fast_magnetosonic_speed(flt_map &prim_vals)
{
  enzo_float cs2 = std::pow(sound_speed(prim_vals),2);
  enzo_float B2 = (prim_vals["bfield_i"]*prim_vals["bfield_i"] +
		   prim_vals["bfield_j"]*prim_vals["bfield_j"] +
		   prim_vals["bfield_k"]*prim_vals["bfield_k"]);
  enzo_float va2 = B2/prim_vals["density"];
  enzo_float cos2 = prim_vals["bfield_i"]*prim_vals["bfield_i"] / B2;
  return std::sqrt(0.5*(va2+cs2+std::sqrt(std::pow(cs2+va2,2) -
					  4.*cs2*va2*cos2)));
}

//----------------------------------------------------------------------

// Applies the pressure_floor to total_energy
void EnzoEOSIdeal::apply_floor_to_energy(Block *block, Grouping &cons_group)
{
  EnzoFieldArrayFactory array_factory(block);
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
