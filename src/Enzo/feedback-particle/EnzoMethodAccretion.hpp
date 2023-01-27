// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretion.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       24 February 2022
/// @brief  Implementation of EnzoMethodAccretion, a base class
///         for accretion methods. These methods compute
///         the accretion rate onto sink particles, and change the properties
///         of the particles accordingly. Gas density is reduced by setting values
///         for the "density_accreted" field. The "accretion_remove_gas" method
///         then subtracts density_accreted from the gas density field

#ifndef ENZO_ENZO_METHOD_ACCRETION
#define ENZO_ENZO_METHOD_ACCRETION

class EnzoMethodAccretion : public Method {

  /// @class   EnzoMethodAccretion
  /// @ingroup Enzo
  /// @brief   [\ref Enzo] Base class for accretion methods, which add
  ///          mass and momentum to sink particles while removing it from the gas
  ///          within an accretion zone around the sink particle.

public:

  // Constructor
  EnzoMethodAccretion(double accretion_radius_cells,
		      double physical_density_threshold_cgs,
		      double max_mass_fraction);

  /// Destructor
  virtual ~EnzoMethodAccretion() throw() {};

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodAccretion);

  /// Charm++ PUP::able migration constructor
  EnzoMethodAccretion (CkMigrateMessage *m)
    : Method (m)
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method
  virtual void compute ( Block * block) throw();

  /// Name
  virtual std::string name () throw()
  { return "accretion"; }

  /// Not sure if this is needed
  virtual std::string particle_type () throw()
  { return "sink";}

  // Compute the maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

  // Update fields representing fluid quantities.
  void update_fields(EnzoBlock * enzo_block) throw();


protected:

  // methods

  // Checks to be performed at initial cycle
  void do_checks_(const Block* block) throw();

  // attributes

  // The accretion radius relative to the cell width
  double accretion_radius_cells_;

  // `physical_density_threshold_cgs_` has a different interpretation
  // depending on the flavor of accretion used. In all
  // cases, it has to be at least as large as the
  // density floor imposed by the VL+CT method.
  //
  // "physical_density_threshold_cgs" accretion: for all cells in
  // accretion zone, if the density in a cell is higher
  // than physical_density_threshold_cgs_, the density is reduced to
  // max(physical_density_threshold_cgs_,(1-max_mass_fraction_)*density_)

  // For other accretion methods, it is the minimum value
  // that the gas density in any cell can take after accretion.
  double physical_density_threshold_cgs_;

  // The maximum fraction of mass that can be accreted from a cell
  // in one timestep
  double max_mass_fraction_;

  // ID for this method's "accumulate refresh"
  int ir_accretion_;

};

#endif // ENZO_ENZO_METHOD_ACCRETION
