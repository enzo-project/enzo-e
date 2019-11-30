// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSourceInternalEnergy.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri December 13 2019
/// @brief    [\ref Enzo] Implementation of Enzo's SourceInternalEnergy class.

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoSourceInternalEnergy::calculate_source
(Block *block, double dt, Grouping &prim_group, Grouping &dUcons_group,
 std::string interface_velocity_name, int dim, EnzoEquationOfState *eos,
 int stale_depth) const throw()
{
  // SANITY CHECKS:
  ASSERT("EnzoSourceInternalEnergy::calculate_source",
	 ("This function should not be called if the dual energy formalism "
	  "is not in use."),
	 eos-> uses_dual_energy_formalism());
  ASSERT("EnzoSourceInternalEnergy::calculate_source",
	 "The EOS can't be barotropic and use the dual energy formalism.",
	 !(eos->is_barotropic()) );

  EnzoBlock * enzo_block = enzo::block(block);
  enzo_float dtdx = dt/enzo_block->CellWidth[dim];

  EnzoFieldArrayFactory array_factory(block,stale_depth);
  EnzoPermutedCoordinates coord(dim);

  // rho and eint are the cell-centered density and specific internal energy
  EFlt3DArray rho = array_factory.from_grouping(prim_group, "density", 0);
  EFlt3DArray eint = array_factory.from_grouping(prim_group,
						 "internal_energy", 0);

  // deint_dens holds the accumulated total change in the internal energy
  // density over the current time-step
  EFlt3DArray deint_dens = array_factory.from_grouping(dUcons_group,
						       "internal_energy", 0);

  // define vl(i) and vr(i) (where i increases in value along dimension dim)
  // as the values of the interface_velocity_field at i-1/2 and i+1/2
  EFlt3DArray vl, vr;
  vl = array_factory.assigned_center_from_name(interface_velocity_name, dim);
  vr = coord.left_edge_offset(vl, 0, 0, 1);

  enzo_float gm1 = eos->get_gamma() - 1.;

  // in the original flux_hll.F the floor was set to a tiny (a small number)
  enzo_float p_floor = eos->get_pressure_floor();

  for (int iz=1; iz<eint.shape(0)-1; iz++) {
    for (int iy=1; iy<eint.shape(1)-1; iy++) {
      for (int ix=1; ix<eint.shape(2)-1; ix++) {
	enzo_float p = gm1*eint(iz,iy,ix)*rho(iz,iy,ix);
	// the following just applies std::max (in a macro-enabled debug mode
	// it will raise errors when the floor is actually needed)
	p = EnzoEquationOfState::apply_floor(p, p_floor);
	deint_dens(iz,iy,ix) += ( dtdx * p * (vl(iz,iy,ix) - vr(iz,iy,ix)) );
      }
    }
  }
}
