#ifndef ENZO_ENZO_INITIAL_ISOLATED_GALAXY_HPP
#define ENZO_ENZO_INITIAL_ISOLATED_GALAXY_HPP

class EnzoInitialIsolatedGalaxy : public Initial {

public: // interface

  /// CHARM++ constructor
  EnzoInitialIsolatedGalaxy(const EnzoConfig * enzo_config) throw();

  // do I need to init things here????
  PUPable_decl(EnzoInitialIsolatedGalaxy); // do i need this??? AE

  /// Charm++ PUP::able migration constructor (AE unsure of this?)
  EnzoInitialIsolatedGalaxy (CkMigrateMessage *m)
    : Initial (m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize block
  virtual void enforce_block
  ( Block * block,
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy * hierarchy ) throw();

  /// Read in particle data (DM and stars)
  void ReadInParticleData(void);

  /// Read in circular velocity table
  void    ReadInVcircData(void);

  /// Intarpolate circular velocity table
  double InterpolateVcircTable(double radius);

  /// Compute cell mass with gaussian interpolation
  double gauss_mass(const double rho_zero,
                    const double xpos,
                    const double ypos,
                    const double zpos,
                    const double dx);

  /// Destructor
  virtual ~EnzoInitialIsolatedGalaxy() throw() {};

private: // attributes

  double center_position_[3];
  double bfield_[3];
  double scale_length_;
  double scale_height_;
  double disk_mass_;
  double gas_fraction_;
  double disk_temperature_;
  double disk_metallicity_;
  double gas_halo_mass_;
  double gas_halo_temperature_;
  double gas_halo_metallicity_;

  double gas_halo_density_;
  double gas_halo_radius_;

  int dual_energy_;
  double uniform_density_;
  double gamma_;
  double mu_;

  const int VCIRC_TABLE_LENGTH = 10000;

  double vcirc_radius[10000];
  double vcirc_velocity[10000];

  double mass_units_;
  double length_units_;
  double density_units_;
  double time_units_;

};

#endif /* ENZO_ENZO_INITIAL_ISOLATED_GALAXY_HPP */
