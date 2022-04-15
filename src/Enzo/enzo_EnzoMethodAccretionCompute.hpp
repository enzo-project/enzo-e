// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretionCompute.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       24 February 2022
/// @brief  Implementation of EnzoMethodAccretionCompute, a base class
///         for "accretion compute" methods. These methods compute
///         the accretion rate onto star particles, and change the properties
///         of the particles accordingly. Gas density is reduced by setting values
///         for the "density_accreted" field. The "accretion_remove_gas" method
///         then subtracts density_accreted from the gas density field

#ifndef ENZO_ENZO_METHOD_ACCRETION_COMPUTE
#define ENZO_ENZO_METHOD_ACCRETION_COMPUTE

class EnzoMethodAccretionCompute : public Method {

  /// @class   EnzoMethodAccretionCompute
  /// @ingroup Enzo
  /// @brief   [\ref Enzo] Base class for AccretionCompute methods, which add
  ///          mass to star particles and set negative values for density_accreted
  ///          field within an accretion zone.

public:

  // Constructor
  EnzoMethodAccretionCompute(double accretion_radius_cells,
			     double density_threshold,
			     double max_mass_fraction,
			     bool conserve_angular_momentum);

  /// Destructor
  virtual ~EnzoMethodAccretionCompute() throw() {};

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodAccretionCompute);

  /// Charm++ PUP::able migration constructor
  EnzoMethodAccretionCompute (CkMigrateMessage *m)
    : Method (m)
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method
  virtual void compute ( Block * block) throw();

  /// Name
  virtual std::string name () throw()
  { return "accretion_compute"; }

  /// Not sure if this is needed
  virtual std::string particle_type () throw()
  { return "star";}

  // Compute the maximum timestep for this method
  virtual double timestep ( Block * block) const throw();


protected:

  // methods

  // Checks to be performed at initial cycle
  void do_checks_() throw();

  // attributes

  // The accretion radius relative to the cell width
  double accretion_radius_cells_;

  // density_threshold_ has a different interpretation
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
  // (with respect to the star particle) is conserved during accretion
  // (See Bleuler and Teyssier 2014, MNRAS 445, 4015–4036 and
  // Krumholz+ 2004, ApJ, 611, 399 for details).
  bool conserve_angular_momentum_;

};

#endif /* EnzoMethodAccretionCompute */
