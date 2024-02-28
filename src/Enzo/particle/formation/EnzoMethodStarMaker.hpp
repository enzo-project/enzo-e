// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMaker.hpp
/// @author     Andrew Emerick (aemerick11@gmail.com)
/// @date
/// @brief

#ifndef ENZO_ENZO_METHOD_STARMAKER
#define ENZO_ENZO_METHOD_STARMAKER

///
/// Game plan. Develop a general star maker class
///   (and an analagous star feedback class?) that
///   will be used to make specific classes for various methods
///   this *may* help things in the future? At least it logically
///   makes sense to do something like this

class EnzoMethodStarMaker : public Method {

  /// @class   EnzoMethodStarMaker
  /// @ingroup Enzo
  /// @btief   [\ref Enzo] Encapsulate Star Maker Routines

public:

  // Create a new StarMaker object
  EnzoMethodStarMaker();

  /// Destructor
  virtual ~EnzoMethodStarMaker() throw() {};

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodStarMaker);

  /// Charm++ PUP::able migration constructor
  EnzoMethodStarMaker (CkMigrateMessage *m)
    : Method (m)
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method
  virtual void compute( Block * block) throw();

  /// Name
  virtual std::string name () throw()
  { return "star_maker"; }

  /// Default particle type for star maker -
  ///   this should be overwritten in child
  ///   classes that may not be making normal stars
  virtual std::string particle_type () throw()
  { return "star";}

//  static void convert_densities_to_fraction(EnzoBlock * enzo_block,
//                                            int direction) throw();

  // Function to rescale fraction fields to new density -
  //   AE: NOTE: This could likely be moved to be a function of
  //             a more general class, not specific to star maker
  static void rescale_densities(EnzoBlock * enzo_block,
                                const int index,
                                const double density_ratio) throw();

protected: // methods

  // Routine functions for checking certain conditions
  //
  int check_number_density_threshold(const double &d);
  int check_overdensity_threshold(const double &rho);
  int check_velocity_divergence(
                enzo_float *vx, enzo_float *vy, enzo_float *vz,
                const int &index, const int &dix, const int &diy,
                const int &diz,
                const double dx, const double dy, const double dz);

  int check_mass(const double &m);

  int check_self_gravitating(
    const double mean_particle_mass, const double density,
    const enzo_float temperature,
    enzo_float *vx, enzo_float *vy, enzo_float *vz,
    const double lunit, const double vunit, const double rhounit, 
    const int &index, const int &dix, const int &diy, const int &diz,
    const double dx, const double dy, const double dz);

  int check_self_gravitating_alt(const double total_energy, const double potential);
  double h2_self_shielding_factor(
    enzo_float *rho, const double metallicity,
    const double dunit, const double lunit,
    const int &index, const int &dix, const int &diy, const int &diz,
    const double dx, const double dy, const double dz);
  int check_jeans_mass(
    const double temperature, const double mean_particle_mass,
    const double density, const double mass,
    const double munit, const double rhounit);
  int check_cooling_time(const double &cooling_time, const double &total_density,
    const double rhounit, const double tunit);
  int check_metallicity(const double &Z);
  int check_temperature(const double &T);


protected: // attributes

  bool use_density_threshold_;
  bool use_velocity_divergence_;
  bool use_self_gravitating_;
  bool use_altAlpha_;
  bool use_h2_self_shielding_;
  bool use_jeans_mass_;
  bool use_overdensity_threshold_;
  bool use_critical_metallicity_;
  bool use_temperature_threshold_;
  bool use_cooling_time_;
  bool use_dynamical_time_;
  double overdensity_threshold_;
  double critical_metallicity_;
  double number_density_threshold_;
  double efficiency_;
  double maximum_star_fraction_;
  double star_particle_min_mass_;
  double star_particle_max_mass_;
  double temperature_threshold_;
  // variables to be passsed here
};

#endif /* EnzoMethodStarMaker */
