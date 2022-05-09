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
		      double density_threshold,
		      double max_mass_fraction,
		      bool conserve_angular_momentum,
		      double angular_momentum_threshold);

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

  // Update fields representing fluid quantities with the "source fields".
  void add_source_fields(EnzoBlock * enzo_block) throw();


protected:

  // methods

  // Checks to be performed at initial cycle
  void do_checks_() throw();

  // attributes

  // The accretion radius relative to the cell width
  double accretion_radius_cells_;

  // `density_threshold_` has a different interpretation
  // depending on the flavor of accretion used. In all
  // cases, it has to be at least as large as the
  // density floor imposed by the VL+CT method.
  //
  // "density_threshold" accretion: for all cells in
  // accretion zone, if the density in a cell is higher
  // than density_threshold_, the density is reduced to
  // max(density_threshold_,(1-max_mass_fraction_)*density_)

  // For other accretion methods, it is the minimum value
  // that the gas density in any cell can take after accretion.
  double density_threshold_;

  // The maximum fraction of mass that can be accreted from a cell
  // in one timestep
  double max_mass_fraction_;

  // If true, angular momentum of the gas in the accretion zone
  // (with respect to the sink particle) is conserved during accretion
  // (See Bleuler and Teyssier 2014, MNRAS 445, 4015–4036 and
  // Krumholz+ 2004, ApJ, 611, 399 for details).
  bool conserve_angular_momentum_;

  // This attribute specifies the angular momentum threshold radius in terms
  // of the minimum cell width at the highest level of refinement.
  // If `conserve_angular_momentum_` is true, then the angular momentum of the gas
  // in a given cell is conserved if and only if its distance from the
  // sink particle is greater than the angular momentum threshold radius.
  // If `conserve_angular_momentum_` is false, this attribute is ignored.
  double ang_mom_threshold_radius_cells_;

  // ID for this method's "accumulate refresh"
  int ir_accretion_;

};

#endif // ENZO_ENZO_METHOD_ACCRETION
