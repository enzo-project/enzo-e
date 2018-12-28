
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

void EnzoEOSIdeal::primitive_from_conservative(Block *block,
					       Grouping &cons_group,
					       Grouping &prim_group)
{
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();
  // get integration limits 
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];

  enzo_float *cons_dens, *px, *py, *pz, *etot, *cons_bx, *cons_by, *cons_bz;
  cons_dens = load_grouping_field_(&field,&cons_group,"density",0);
  px = load_grouping_field_(&field,&cons_group,"momentum",0);
  py = load_grouping_field_(&field,&cons_group,"momentum",1);
  pz = load_grouping_field_(&field,&cons_group,"momentum",2);
  etot = load_grouping_field_(&field,&cons_group,"total_energy",0);
  cons_bx = load_grouping_field_(&field,&cons_group,"bfield",0);
  cons_by = load_grouping_field_(&field,&cons_group,"bfield",1);
  cons_bz = load_grouping_field_(&field,&cons_group,"bfield",2);

  enzo_float *prim_dens, *vx, *vy, *vz, *pressure, *prim_bx, *prim_by, *prim_bz;
  prim_dens = load_grouping_field_(&field,&prim_group,"density",0);
  vx = load_grouping_field_(&field,&prim_group,"velocity",0);
  vy = load_grouping_field_(&field,&prim_group,"velocity",1);
  vz = load_grouping_field_(&field,&prim_group,"velocity",2);
  pressure = load_grouping_field_(&field,&prim_group,"pressure",0);
  prim_bx = load_grouping_field_(&field,&prim_group,"bfield",0);
  prim_by = load_grouping_field_(&field,&prim_group,"bfield",1);
  prim_bz = load_grouping_field_(&field,&prim_group,"bfield",2);

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
	int i = ix + mx*(iy + my*iz);

	enzo_float density = cons_dens[i];
	prim_dens[i] = density;
	vx[i] = px[i]/density;
	vy[i] = py[i]/density;
	vz[i] = pz[i]/density;
	enzo_float v2 = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	pressure[i] = (etot[i] - 0.5*v2*density)*(gamma_-1.);
	
	prim_bx[i] = cons_bx[i];
	prim_by[i] = cons_by[i];
	prim_bz[i] = cons_bz[i];
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::conservative_from_primitive(flt_map &prim, flt_map &cons){
  cons["density"] = prim["density"];
  cons["momentum_i"] = prim["density"]*prim["velocity_i"];
  cons["momentum_j"] = prim["density"]*prim["velocity_j"];
  cons["momentum_k"] = prim["density"]*prim["velocity_k"];
  enzo_float v2 = (prim["velocity_i"]*prim["velocity_i"] +
		   prim["velocity_j"]*prim["velocity_j"] +
		   prim["velocity_k"]*prim["velocity_k"]);
  enzo_float e_int = prim["pressure"]/(gamma_-1.);
  cons["total_energy"] = e_int + 0.5*v2*prim["density"];
  cons["bfield_i"] = prim["bfield_i"];
  cons["bfield_j"] = prim["bfield_j"];
  cons["bfield_k"] = prim["bfield_k"];
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
