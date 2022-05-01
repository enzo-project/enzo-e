// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoPhysicsFluidProps.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-04-05
/// @brief    [\ref Enzo] Implementation of the EnzoPhysicsFluidProps class

#include "enzo.hpp"

//----------------------------------------------------------------------

// this is to be replaced in the future!
namespace{

  bool has_method_(const std::string& name){
    const EnzoConfig * config = enzo::config();
    const std::vector<std::string>& mlist = config->method_list;
    return std::find(mlist.begin(), mlist.end(), name) != mlist.end();
  }

  bool uses_ppm_() { return has_method_("ppm"); }
  bool uses_vlct_() { return has_method_("mhd_vlct"); }

  EnzoDualEnergyConfig get_de_config_(){
    const EnzoConfig * config = enzo::config();

    // these use the legacy configuration options!
    if (uses_ppm_() & config->ppm_dual_energy){
      ASSERT("get_deconfig_",
             "a simulation can't use ppm and vlct solvers at once",
             !uses_vlct_());
      return EnzoDualEnergyConfig::build_bryan95_formulation
        (config->ppm_dual_energy_eta_1,
         config->ppm_dual_energy_eta_2);
    } else if (uses_vlct_() & config->method_vlct_dual_energy){
      return EnzoDualEnergyConfig::build_modern_formulation
        (config->method_vlct_dual_energy_eta);
    }
    return EnzoDualEnergyConfig::build_disabled();
  }

  EnzoFluidFloorConfig get_fluid_floor_config_(){
    const EnzoConfig * config = enzo::config();

    double density_floor = -1;
    double pressure_floor = -1;
    double temperature_floor = -1;

    // these use the legacy configuration options!
    if (uses_ppm_()){
      ASSERT("get_fluid_floor_config_",
             "a simulation can't use ppm and vlct solvers at once",
             !uses_vlct_());
      density_floor = config->ppm_density_floor;
      pressure_floor = config->ppm_pressure_floor;
      temperature_floor = config->ppm_temperature_floor;
    } else if (uses_vlct_()){
      density_floor = config->method_vlct_density_floor;
      pressure_floor = config->method_vlct_pressure_floor;
    }
    return {density_floor, pressure_floor, temperature_floor};
  }

}

EnzoPhysicsFluidProps::EnzoPhysicsFluidProps() noexcept
  : Physics(),
    de_config_(get_de_config_()),
    fluid_floor_config_(get_fluid_floor_config_()),
    gamma_(0.0),
    mol_weight_(0.0)
{
  const EnzoConfig * config = enzo::config();
  gamma_ = config->field_gamma;
  mol_weight_ = config->ppm_mol_weight;
}
