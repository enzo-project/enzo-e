// See LICENSE_CELLO file for license and copyright information
/// @file     enzo_EnzoInitialInclinedWave.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri May 24 2019
/// @brief    [\ref Enzo] Initialization routine for spherical cloud in a wind

#ifndef ENZO_ENZO_INITIAL_CLOUD_HPP
#define ENZO_ENZO_INITIAL_CLOUD_HPP

// tolerance to be used when checking the equivalence of 2 floating point
// numbers in EnzoInitialCloud. This is the same as ETA_TOLERANCE for doubles.
#define INIT_CLOUD_TOLERANCE 1.0e-10

// represents a spherical region
class SphereRegion;

class EnzoInitialCloud : public Initial {
  /// @class    EnzoInitialCloud
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer a cloud crushing problem
  /// The wind is assumed to blow along the x-directions

public: // interface

  /// Constructor
  EnzoInitialCloud(int cycle, double time, int subsample_n,
		   double cloud_radius, double cloud_center_x,
		   double cloud_center_y, double cloud_center_z,
		   double density_cloud, double density_wind,
		   double etot_wind, double eint_wind,
		   double velocity_wind, double metal_mass_frac,
		   double perturb_stddev, double truncate_dev,
		   unsigned int perturb_seed)
    : Initial(cycle,time),
      subsample_n_(subsample_n),
      cloud_radius_(cloud_radius),
      cloud_center_x_(cloud_center_x),
      cloud_center_y_(cloud_center_y),
      cloud_center_z_(cloud_center_z),
      density_cloud_(density_cloud),
      density_wind_(density_wind),
      etot_wind_(etot_wind),
      eint_wind_(eint_wind),
      velocity_wind_(velocity_wind),
      metal_mass_frac_(metal_mass_frac),
      perturb_stddev_(perturb_stddev),
      truncate_dev_(truncate_dev),
      perturb_seed_(perturb_seed)
  {
    ASSERT("EnzoInitialCloud", "subsample_n must be >=0", subsample_n>=0);
    ASSERT("EnzoInitialCloud", "cloud_radius must be positive",
	   cloud_radius>0);
    ASSERT("EnzoInitialCloud", "density_wind must be positive",
           density_wind_>0);
    ASSERT("EnzoInitialCloud", "density_cloud must be positive",
           density_cloud_>0);
    ASSERT("EnzoInitialCloud", "etot_wind must exceed wind kinetic energy",
	   etot_wind>0.5*velocity_wind_*velocity_wind_);
    ASSERT("EnzoInitialCloud", "eint_wind must be zero or positive.",
	   eint_wind>=0.);
    ASSERT("EnzoInitialCloud", "metal_mass_frac must be in [0,1]",
	   metal_mass_frac >=0 && metal_mass_frac <=1);
    ASSERT("EnzoInitialCloud", "perturb_stddev must be >= 0",
	   perturb_stddev >=0);
    ASSERT("EnzoInitialCloud", "truncate_dev must be >= 0", truncate_dev >=0);
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialCloud);

  /// CHARM++ migration constructor
  EnzoInitialCloud(CkMigrateMessage *m)
    : Initial(m),
      subsample_n_(0),
      cloud_radius_(0.),
      cloud_center_x_(0.),
      cloud_center_y_(0.),
      cloud_center_z_(0.),
      density_cloud_(0.),
      density_wind_(0.),
      etot_wind_(0.),
      velocity_wind_(0.),
      metal_mass_frac_(0.)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: update whenever attributes change

    TRACEPUP;

    Initial::pup(p);
    p | subsample_n_;
    p | cloud_radius_;
    p | cloud_center_x_;
    p | cloud_center_y_;
    p | cloud_center_z_;
    p | density_cloud_;
    p | density_wind_;
    p | etot_wind_;
    p | velocity_wind_;
    p | metal_mass_frac_;
    p | perturb_stddev_;
    p | truncate_dev_;
    p | perturb_seed_;
  }

public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// number of subsampled cells per cell side is (int)pow(2,subsample_n_);
  int subsample_n_;
  /// Cloud radius in code units
  double cloud_radius_;
  /// Cloud center
  double cloud_center_x_;
  double cloud_center_y_;
  double cloud_center_z_;

  /// cloud and wind densities
  double density_cloud_;
  double density_wind_;

  /// specific total energy of the wind
  double etot_wind_;

  /// specific internal energy of the wind
  double eint_wind_;

  /// velocity of the wind
  double velocity_wind_;

  double metal_mass_frac_;

  /// The stddev for normal distribution used to perturb cloud density
  double perturb_stddev_;

  /// Number of stddevs where normal distribution is truncated
  double truncate_dev_;

  /// The random seed used to seed the density perturbations
  unsigned int perturb_seed_;
};

#endif //ENZO_ENZO_INITIAL_CLOUD_HPP
