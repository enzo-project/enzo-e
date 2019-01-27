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

  // Compute the maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

//  static void convert_densities_to_fraction(EnzoBlock * enzo_block,
//                                            int direction) throw();

  void rescale_densities(EnzoBlock * enzo_block,
                         const int index,
                         const double density_ratio) throw();

protected: // methods

  // Routine functions for checking certain conditions
  int check_number_density_threshold(const double &d);
  int check_velocity_divergence(
                double *vx, double *vy, double *vz,
                const int &index, const int &dix, const int &diy,
                const int &diz);
  int check_minimum_mass(const double &m);

protected: // attributes

  bool use_density_threshold_;
  bool use_velocity_divergence_;
  bool use_dynamical_time_;
  double number_density_threshold_;
  double efficiency_;
  double maximum_star_fraction_;
  double star_particle_mass_;
  // variables to be passsed here
};

#endif /* EnzoMethodStarMaker */
