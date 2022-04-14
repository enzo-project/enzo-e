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
			     double density_threshold);

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
  // the value of density_threshold_, and the removed
  // mass and momentum is added to the star particle
  //
  // "bondi_hoyle" accretion: the value of
  // density_threshold_ sets a lower limit for the gas
  // density field after accretion, and so in effect puts
  // a limit on how much gas a star particle can accrete.
  double density_threshold_;

};

#endif /* EnzoMethodAccretionCompute */
