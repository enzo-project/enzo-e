// This class is composed of several component classses
// These classes include the following:
//
//    EnzoEquationOfState:       Equation of State for Gas
//    EnzoReconstructor:         Reconstructs primitive variables
//    EnzoRiemann:               Solves the Riemann Problem
//    EnzoConstrainedTransport:  Performs constrained transport
//
// Some notes on implementation
//    - this Method tracks total energy density, (referred to as total_energy)
//
//    Grouping Objects
//    ----------------
//    Implementation relies upon passing around instances of the Grouping
//    object between different object methods. Originally attempted to pass
//    around vectors of field ids, but that doesn't offer nearly as much
//    flexibility. In particular:
//        - It would make it more difficult to add new fields without altering
//          the interface (e.g. optionally using a dual energy formalism,
//          species, colors, Cosmic Ray energy density and flux density)
//        - Relatedly, by sharing common group names between different grouping
//          objects it is easier to load in relevant fields representing
//          different physical quantities (e.g. Groupings of fluxes and
//          conserved variables each have a "momentum" group)
//
//    Notes about current groupings:
//        - Currently track 14 different groupings:
//            1. conserved quantities
//            2. interface b-fields (longitudinal B-fields centered at
//               interface between cells)
//            3. primitive quantities
//            4. reconstructed left primitive fields
//            5. reconstructed right primitive fields
//            6. fluxes in the x-direction
//            7. fluxes in the y-direction
//            8. fluxes in the z-direction
//            9. face centered fields along x,y,z direction that identify which
//               direction is upwind. (known as weight_group)
//           10. edge centered electric fields
//           11. temp conserved quantities (to store values at the half
//               timestep)
//           12. temp interface b-fields (to store values at the half timestep)
//           13. reconstructed left conserved fields
//           14. reconstructed right conserved fields
//        - several groups only contain 1 field (e.g. the "density" and
//          "total_energy" groups). In effect, they serve as alias names for
//          the fields
//        - This implementation relies on the fact that the fields in a
//          Grouping are sorted alphabetically (in this way the fields in the
//          "momentum" group are always ordered "momentum_x","momentum_y",
//          "momentum_z"
//        - The conserved quantities grouping, x/y/z flux groupings, and
//          temp conserved quantities grouping each always contain groups
//          called "density", "momentum", "total_energy", "bfield"
//        - All fields within the primitive quantites, conserved quanties, and
//          temp conserved quantities are all cell-centered.
//        - All face-centered fields with the exception of (temporary)
//          interface b-fields only include space for the interior of the grid.
//          The (temporary) interface bfields all include space for values on
//          the exterior of the grid.
//        - All reconstructed fields are not technically registered as
//          cell-centered fields. They are reused to store reconstructed values
//          along different axes. Consequently, they are formally registerred as
//          cell-centered fields (to guaruntee that they have enough space).
//        - reconstructed left and right conserved fields are entirely
//          non-essential. They only exist to allow for symmetry between the
//          EquationOfState interface methods: primitive_from_conservative and
//          conservative_from_primitive. (They both operate on groupings of
//          fields). There are three main alternatives to using these groupings:
//             1. Adapt Riemann Solver to use modified flux formulas that don't
//                explicitly require knowledge of the interface flux values
//                (This may couple the RiemannSolver to the gas physics)
//             2. Come up with an elegant and easily extendible way for the
//                EquationOfState to convert primitives to conserved quantities
//                at a single location.
//             3. Completely couple RiemannSolvers to the gas physics and have
//                write out the conversion within the function.
//        - The current implementation of Groupings will not scale once colors
//          and species are implemented.
//
//    The number of tracked Groupings could be reduced drastically if we could
//    track histories of a subset of fields temporarily, and if we could easily
//    search for fields by specifying multiple group tags.
//
//    Mapping
//    -------
//    Current implementation relies on EnzoEquationOfState converting the
//    entire grid of conserved quantities to primitives. Within the Riemann
//    Solver and the instance methods of EnzoMethodVlct (to compute timestep and
//    to add add flux divergence), subclasses of EnzoEquationOfState are used
//    to compute properties of individual cells. These values are passed via a
//    map. We alias the map of values at a single location as a flt_map. The
//     current map implementation is an
//    unordered_map
//        - There is an inconsistency between the maps and groupings. Groupings
//          for vector quantites (i.e. velocity/momentum/bfields) only have a
//          single name, which then corresponds to collection of components. In
//          contrast, the map has an individual key for each component.
//          (e.g. velocity_i, velocity_j, velocity_k)
//        - Can probably come up with a more elegant solution for Mapping and
//          Grouping (the obvious choice is to do away with Mapping completely
//          and just apply operations on Groupings)
//    The Mapping is hard to debug. If you attempt to access a field that does
//    not already exist, a default value will be returned. (0.0 for a float)
//  Issue: Concerned about the units of the magnetic field (especially during
//  calculations of the Alfven Velocity)

#ifndef ENZO_ENZO_METHOD_VLCT_HPP
#define ENZO_ENZO_METHOD_VLCT_HPP

class EnzoEquationOfState;
class EnzoReconstructor;
class EnzoRiemann;
class EnzoConstrainedTransport;

// Frequently used Helper class:
class EnzoPermutedCoordinates;

// for brevity, and convenience (in case the choice of map is altered)
// we define an alias for our choice of map
typedef std::unordered_map<std::string,enzo_float> flt_map;

class EnzoMethodVlct : public Method {

  /// @class    EnzoMethodVlct
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate VL + CT MHD method

public: // interface

  /// Create a new EnzoMethodVlct object
  EnzoMethodVlct(double gamma);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodVlct);

  /// Charm++ PUP::able migration constructor
  EnzoMethodVlct (CkMigrateMessage *m)
    : Method (m),
      eos_(NULL),
      half_dt_recon_(NULL),
      full_dt_recon_(NULL),
      riemann_solver_(NULL)
  {
    // We want to pack/unpack Groupings in the future.
    conserved_group_ = new Grouping;
    bfieldi_group_ = new Grouping;
    primitive_group_ = new Grouping;
    setup_groups_();
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Delete EnzoMethodVlct object
  ~EnzoMethodVlct();

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "vlct"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();


  static std::vector<std::string> cons_group_names;
  static std::vector<std::string> prim_group_names;
protected: // methods

  // Prepare the group
  void setup_groups_();

  // not sure if I will pass field_ids and blocks or arrays
  // not sure if this should be static
  void compute_flux_(Block *block, int dim, Grouping &cur_cons_group,
		     Grouping &cur_bfieldi_group, Grouping &priml_group,
		     Grouping &primr_group, Grouping &flux_group,
		     Grouping &consl_group, Grouping &consr_group, 
		     Grouping &weight_group, EnzoReconstructor &reconstructor);

  // compute the Electric fields using the fluxes and cell-centered
  // primitives
  void compute_efields_(Block *block, Grouping &xflux_group,
			Grouping &yflux_group, Grouping &zflux_group,
			std::string center_efield_name, Grouping &efield_group,
			Grouping &weight_group,
			EnzoConstrainedTransport &ct);

  // adds flux divergence to the fields listed in conserved_group_ and stores
  // the results in out_cons_group (this can be the same as conserved_group_)
  void update_quantities_(Block *block, Grouping &xflux_group,
			  Grouping &yflux_group, Grouping &zflux_group,
			  Grouping &out_cons_group, double dt);

  // allocate the temporary fields needed for scratch space and store their ids
  // efield_ids refer to efields centered on the edges of cells
  void allocate_temp_fields_(Block *block, Grouping &priml_group,
			     Grouping &primr_group, Grouping &xflux_group,
			     Grouping &yflux_group, Grouping &zflux_group,
			     Grouping &efield_group,
			     std::string &center_efield_name,
			     Grouping &weight_group,
			     Grouping &temp_conserved_group,
			     Grouping &temp_bfieldi_group,
			     Grouping &consl_group, Grouping &consr_group);

  // deallocate the temporary fields used for scratch space
  void deallocate_temp_fields_(Block *block, Grouping &priml_group,
			       Grouping &primr_group, Grouping &xflux_group,
			       Grouping &yflux_group, Grouping &zflux_group,
			       Grouping &efield_group,
			       std::string center_efield_name,
			       Grouping &weight_group,
			       Grouping &temp_conserved_group,
			       Grouping &temp_bfieldi_group,
			       Grouping &consl_group, Grouping &consr_group);

protected: // attributes

  // In a lot of senses, these will serve as aliases
  Grouping *conserved_group_;
  Grouping *bfieldi_group_;
  Grouping *primitive_group_;

  EnzoEquationOfState *eos_;
  EnzoReconstructor *half_dt_recon_;
  EnzoReconstructor *full_dt_recon_;
  EnzoRiemann *riemann_solver_;

};

#endif /* ENZO_ENZO_METHOD_VLCT_HPP */
