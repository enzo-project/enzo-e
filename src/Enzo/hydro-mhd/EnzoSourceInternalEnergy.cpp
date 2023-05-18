// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSourceInternalEnergy.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri December 13 2019
/// @brief    [\ref Enzo] Implementation of Enzo's SourceInternalEnergy class.

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoSourceInternalEnergy::calculate_source
(const int dim, const double dt, const enzo_float cell_width,
 const EnzoEFltArrayMap &prim_map, EnzoEFltArrayMap &dUcons_map,
 const CelloArray<const enzo_float,3> &interface_velocity,
 const EnzoEquationOfState *eos,
 const int stale_depth) const throw()
{
  // SANITY CHECKS:
  ASSERT("EnzoSourceInternalEnergy::calculate_source",
	 ("This function should not be called if the dual energy formalism "
	  "is not in use."),
	 eos-> uses_dual_energy_formalism());
  ASSERT("EnzoSourceInternalEnergy::calculate_source",
	 "The EOS can't be barotropic and use the dual energy formalism.",
	 !(eos->is_barotropic()) );

  const enzo_float dtdx = dt/cell_width;

  EnzoPermutedCoordinates coord(dim);

  // Since we don't know the ith component of the velocity on the external
  // faces of the outermost layer of cells, we frame the problem as follows:
  //  deint_dens(k,j,i+1) -=
  //     dtdx_i * pressure(k,j,i+1) * (vel_i(k,j,i+3/2) - vel_i(k,j,i+1/2))
  //
  // Note  - deint_dens holds the accumulated change in the internal energy
  //         density over the current time-step
  //       - we compute pressure from the cell-centered density ane specific
  //         internal from the start of the time-step (adopting convention from
  //         flux_hll.F - we also could just use the precomputed field)

  CSlice full_ax(nullptr, nullptr);

  // load cell-centered quantities.
  // define: pressure_center(k,j,i)    ->  pressure(k,j,i+1)
  //         deint_dens_center(k,j,i)  ->  deint_dens(k,j,i+1)
  const CelloArray<const enzo_float,3> pressure = prim_map.get
    ("pressure", stale_depth);
  const CelloArray<const enzo_float,3> pressure_center = coord.get_subarray
    (pressure, full_ax, full_ax, CSlice(1, -1));

  const CelloArray<enzo_float,3> deint_dens = dUcons_map.get
    ("internal_energy", stale_depth);
  const CelloArray<enzo_float,3> deint_dens_center = coord.get_subarray
    (deint_dens, full_ax, full_ax, CSlice(1, -1));

  // remove staled zone from interface_velocity
  CSlice stale_slice(stale_depth,-stale_depth);
  const CelloArray<const enzo_float,3> vel = interface_velocity.subarray
    (stale_slice, stale_slice, stale_slice);
  // load the face-centered values of the ith velocity component
  // define: vl(k,j,i)                 ->  vel_i(k,j,i+1/2)
  //         vr(k,j,i)                 ->  vel_i(k,j,i+3/2)
  const CelloArray<const enzo_float,3> vl = coord.get_subarray
    (vel, full_ax, full_ax, CSlice(0, -1));
  const CelloArray<const enzo_float,3> vr = coord.get_subarray
    (vel, full_ax, full_ax, CSlice(1, nullptr));

  // in the original flux_hll.F the floor was set to tiny (a small number)
  enzo_float p_floor = eos->get_pressure_floor();

  enzo_float gm1 = eos->get_gamma() - 1.;

  for (int iz=0; iz<pressure_center.shape(0); iz++) {
    for (int iy=0; iy<pressure_center.shape(1); iy++) {
      for (int ix=0; ix<pressure_center.shape(2); ix++) {
	enzo_float p = pressure_center(iz,iy,ix);
	// the following just applies std::max (in a macro-enabled debug mode
	// it will raise errors when the floor is actually needed)
        // It's probably redudant to apply the floor here
	p = EnzoEquationOfState::apply_floor(p, p_floor);
	deint_dens_center(iz,iy,ix) -= dtdx*p*(vr(iz,iy,ix) - vl(iz,iy,ix));
      }
    }
  }
}
