// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoConfig.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-02
/// @brief    [\ref Parameters] Declaration of the EnzoConfig class
///

#ifndef PARAMETERS_ENZO_CONFIG_HPP
#define PARAMETERS_ENZO_CONFIG_HPP

#define MAX_FIELDS      30
#define MAX_FILE_GROUPS 10

class Parameters;

class EnzoConfig : public Config {

  /// @class    EnzoConfig
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Declaration of Enzo configuration class

public: // interface

  /// Constructor
  EnzoConfig() throw();

  /// Destructor
  ~EnzoConfig() throw();

  /// Copy constructor
  EnzoConfig(const EnzoConfig & config) throw();

  /// Assignment operator
  EnzoConfig & operator= (const EnzoConfig & config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoConfig);

  /// CHARM++ migration constructor
  EnzoConfig(CkMigrateMessage *m) : Config (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Read values from the Parameters object
  void read (Parameters * parameters) throw();

public: // attributes

  // NOTE: change pup() function whenever attributes change

  // EnzoMethodPpm

  double                     ppm_density_floor;
  bool                       ppm_diffusion;
  bool                       ppm_dual_energy;
  double                     ppm_dual_energy_eta_1;
  double                     ppm_dual_energy_eta_2;
  int                        ppm_flattening;
  int                        ppm_minimum_pressure_support_parameter;
  double                     ppm_number_density_floor;
  double                     ppm_pressure_floor;
  bool                       ppm_pressure_free;
  bool                       ppm_steepening;
  float                      ppm_temperature_floor;
  bool                       ppm_use_minimum_pressure_support;

  double                     field_gamma;

  // Cosmology (NOT ACCESSED)
  bool                       physics_cosmology;
  double                     physics_cosmology_comoving_box_size;
  double                     physics_cosmology_hubble_constant_now;
  double                     physics_cosmology_initial_redshift;
  double                     physics_cosmology_max_expansion_rate;
  double                     physics_cosmology_omega_lamda_now;
  double                     physics_cosmology_omega_matter_now;

  // EnzoInitialSedovArray[23]
  int                        initial_sedov_array[3];
  double                     initial_sedov_radius_relative;
  double                     initial_sedov_pressure_in;
  double                     initial_sedov_pressure_out;
  double                     initial_sedov_density;

  // EnzoProlong
  std::string                interpolation_method;

  // EnzoMethodHeat
  double                     method_heat_alpha;

  // EnzoMethodNull
  double                     method_null_dt;

  // EnzoMethodGrackle

  // Units

  code_units method_grackle_units;
  double method_grackle_gamma;
  bool   method_grackle_with_radiative_cooling;
  bool   method_grackle_primordial_chemistry;
  bool   method_grackle_metal_cooling;
  bool   method_grackle_h2_formation_on_dust;

  bool   method_grackle_cmb_temperature_floor;
  std::string method_grackle_data_file;

  int    method_grackle_three_body_rate;
  int    method_grackle_cie_cooling;
  int    method_grackle_h2_optical_depth_approximation;

  int    method_grackle_photoelectric_heating;
  double method_grackle_photoelectric_heating_rate;

  int    method_grackle_UVbackground;

  double method_grackle_UVbackground_redshift_on;
  double method_grackle_UVbackground_redshift_off;
  double method_grackle_UVbackground_redshift_fullon;
  double method_grackle_UVbackground_redshift_drop;

  int    method_grackle_Compton_xray_heating;

  double method_grackle_LWbackground_intensity;
  int    method_grackle_LWbackground_sawtooth_suppression;

  double method_grackle_HydrogenFractionByMass;
  double method_grackle_DeuteriumToHydrogenRatio;
  double method_grackle_SolarMetalFractionByMass;
  int    method_grackle_NumberOfTemperatureBins;
  int    method_grackle_ih2co;
  int    method_grackle_ipiht;
  double method_grackle_TemperatureStart;
  double method_grackle_TemperatureEnd;
  int    method_grackle_comp_xray;
  int    method_grackle_temp_xray;
  int    method_grackle_CaseBRecombination;
  int    method_grackle_NumberOfDustTemperatureBins;
  double method_grackle_DustTemperatureStart;
  double method_grackle_DustTemperatureEnd;

  double method_grackle_cloudy_electron_fraction_factor;
  
};

#endif /* PARAMETERS_ENZO_CONFIG_HPP */

