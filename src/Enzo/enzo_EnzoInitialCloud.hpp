// See LICENSE_CELLO file for license and copyright information
/// @file     enzo_EnzoInitialInclinedWave.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri May 24 2019
/// @brief    [\ref Enzo] Initialization routine for spherical cloud in a wind

#ifndef ENZO_ENZO_INITIAL_CLOUD_HPP
#define ENZO_ENZO_INITIAL_CLOUD_HPP

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
		   double pressure, double velocity_wind)
    : Initial(cycle,time),
      subsample_n_(subsample_n),
      cloud_radius_(cloud_radius),
      cloud_center_x_(cloud_center_x),
      cloud_center_y_(cloud_center_y),
      cloud_center_z_(cloud_center_z),
      density_cloud_(density_cloud),
      density_wind_(density_wind),
      pressure_(pressure),
      velocity_wind_(velocity_wind)
  {
    ASSERT("EnzoInitialCloud", "subsample_n must be >=0", subsample_n>=0);
    ASSERT("EnzoInitialCloud", "cloud_radius must be positive",
	   cloud_radius>0);
    ASSERT("EnzoInitialCloud", "density_wind must be positive",
           density_wind_>0);
    ASSERT("EnzoInitialCloud", "density_cloud must be positive",
           density_cloud_>0);
    ASSERT("EnzoInitialCloud", "pressure must be positive",
           pressure>0);
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
      pressure_(0.),
      velocity_wind_(0.)
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
    p | pressure_;
    p | velocity_wind_;
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

  /// This should eventually be replaced with total_energy_wind_ and
  /// total_energy_cloud_ (specific energies)
  double pressure_;

  /// velocity of the wind
  double velocity_wind_;
};

#endif //ENZO_ENZO_INITIAL_CLOUD_HPP
