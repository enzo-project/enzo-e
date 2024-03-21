// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProblem.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-03
/// @brief    Implementation of EnzoProblem class
///
///

#include "enzo.hpp"
#include "Enzo/hydro-mhd/hydro-mhd.hpp" // EnzoMethodMHDVlct, EnzoMethodPpm,
                                        // EnzoMethodPpml

#include "Enzo/assorted/assorted.hpp" // misc. Method classes
#include "Enzo/gravity/gravity.hpp" // EnzoMethodGravity
                                    // EnzoMethodBackgroundAcceleration
                                    // EnzoComputeAcceleration
                                    // EnzoSolver* EnzoMatrix*
#include "Enzo/initial/initial.hpp" // lots of initializers
#include "Enzo/io/io.hpp" // EnzoMethodCheck, EnzoInitial{Hdf5,Music}
#include "Enzo/mesh/mesh.hpp" // EnzoProlong, EnzoRefine*, EnzoRestrict*
#include "Enzo/particle/particle.hpp"
#include "Enzo/tests/tests.hpp" // EnzoInitial*Test
#include "Enzo/utils/utils.hpp" // EnzoComputeCicInterp

//----------------------------------------------------------------------

EnzoProblem::EnzoProblem() throw ()
  : Problem()
{
}

//----------------------------------------------------------------------

EnzoProblem::~EnzoProblem() throw ()
{
}

//----------------------------------------------------------------------

void EnzoProblem::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Problem::pup(p);
}

//======================================================================

Boundary * EnzoProblem::create_boundary_
(std::string type,
 int index,
 Config * config,
 Parameters * parameters
 ) throw ()
/// @param type   Type of boundary condition to use
/// @param index  Index of the Boundary object to use
/// @param config Configuration parameters object
/// @param parameters Parameters object (for evaluating expressions)
{

  std::shared_ptr<Mask> mask = nullptr;
  if (config->boundary_mask[index]) {
    std::string param_str = "Boundary:" + config->boundary_list[index] + ":mask";
    Param * param = parameters->param(param_str);
    mask = Mask::create(param,parameters);
  }
  axis_enum axis = (axis_enum) config->boundary_axis[index];
  face_enum face = (face_enum) config->boundary_face[index];

  Boundary * boundary = 0;

  if (       type == "reflecting") {
    boundary = new EnzoBoundary (axis,face,mask,boundary_type_reflecting);
  } else if (type == "outflow") {
    boundary = new EnzoBoundary (axis,face,mask,boundary_type_outflow);
  } else {
    boundary = Problem::create_boundary_(type,index,config,parameters);
  }

  return boundary;
}

//----------------------------------------------------------------------

Initial * EnzoProblem::create_initial_
(
 std::string  type,
 int index,
 Config * config,
 Parameters * parameters
 ) throw ()
{

  //--------------------------------------------------
  // parameter: Initial : cycle
  // parameter: Initial : time
  //--------------------------------------------------

  Initial * initial = 0;

  // move creation of p_accessor up the call stack?
  const std::string root_path =
    ("Initial:" + parameters->list_value_string(index, "Initial:list"));
  ParameterGroup p_group(*parameters, root_path);

  int cycle   = config->initial_cycle;
  double time = config->initial_time;

  const EnzoConfig * enzo_config = enzo::config();

  if (type == "hdf5") {

    initial = new EnzoInitialHdf5
      (cycle,time,
       enzo_config->initial_hdf5_max_level,
       enzo_config->initial_hdf5_format,
       enzo_config->initial_hdf5_blocking,
       enzo_config->initial_hdf5_monitor_iter,
       enzo_config->initial_hdf5_field_files,
       enzo_config->initial_hdf5_field_datasets,
       enzo_config->initial_hdf5_field_coords,
       enzo_config->initial_hdf5_field_names,
       enzo_config->initial_hdf5_particle_files,
       enzo_config->initial_hdf5_particle_datasets,
       enzo_config->initial_hdf5_particle_coords,
       enzo_config->initial_hdf5_particle_types,
       enzo_config->initial_hdf5_particle_attributes);

  } else if (type == "music") {

    initial = new EnzoInitialMusic
      (cycle,time,enzo_config,config->mesh_max_initial_level);

  } else if (type == "implosion_2d") {

    initial = new EnzoInitialImplosion2(cycle,time);

  } else if (type == "sedov_array_2d") {

    initial = new EnzoInitialSedovArray2(enzo_config);

  } else if (type == "sedov_array_3d") {

    initial = new EnzoInitialSedovArray3(enzo_config);

  } else if (type == "sedov_random") {

    initial = new EnzoInitialSedovRandom(enzo_config);

  } else if (type == "sedov") {

    const int rank = enzo_config->initial_sedov_rank;

    ASSERT1 ("EnzoConfig::read()",
	     "Parameter 'Initial:sedov:rank' is %d, but must be set to 2 or 3",
	     rank,  (rank == 2 || rank == 3) );

    if (rank == 2) initial = new EnzoInitialSedovArray2(enzo_config);
    if (rank == 3) initial = new EnzoInitialSedovArray3(enzo_config);

#ifdef CONFIG_USE_GRACKLE
  } else if (type == "grackle_test") {
    initial = new EnzoInitialGrackleTest(enzo_config);
#endif /* CONFIG_USE_GRACKLE */
  } else if (type == "feedback_test") {
    initial = new EnzoInitialFeedbackTest(enzo_config);
  } else if (type == "vlct_bfield") {
    initial = new EnzoInitialBCenter(parameters, cycle, time,
				     enzo_config->initial_bcenter_update_etot);
  } else if (type == "cloud") {
    initial = new EnzoInitialCloud
      (cycle,time,
       enzo_config->initial_cloud_subsample_n,
       enzo_config->initial_cloud_radius,
       enzo_config->initial_cloud_center_x,
       enzo_config->initial_cloud_center_y,
       enzo_config->initial_cloud_center_z,
       enzo_config->initial_cloud_density_cloud,
       enzo_config->initial_cloud_density_wind,
       enzo_config->initial_cloud_etot_wind,
       enzo_config->initial_cloud_eint_wind,
       enzo_config->initial_cloud_velocity_wind,
       enzo_config->initial_cloud_metal_mass_frac,
       enzo_config->initial_cloud_initialize_uniform_bfield,
       enzo_config->initial_cloud_uniform_bfield,
       enzo_config->initial_cloud_perturb_Nwaves,
       enzo_config->initial_cloud_perturb_amplitude,
       enzo_config->initial_cloud_perturb_min_wavelength,
       enzo_config->initial_cloud_perturb_max_wavelength,
       enzo_config->initial_cloud_perturb_seed);
  } else if (type == "collapse") {
    initial = new EnzoInitialCollapse
      (cycle,time,
       enzo_config->initial_collapse_rank,
       enzo_config->initial_collapse_array,
       enzo_config->initial_collapse_radius_relative,
       enzo_config->initial_collapse_particle_ratio,
       enzo_config->initial_collapse_mass,
       enzo_config->initial_collapse_temperature);
  } else if (type == "cosmology") {
    initial = new EnzoInitialCosmology
      (cycle,time,
       enzo::fluid_props()->gamma(),
       enzo_config->initial_cosmology_temperature
       );
  } else if (type == "inclined_wave") {
    initial = new EnzoInitialInclinedWave(cycle, time, p_group);
  } else if (type == "turbulence") {
    initial = new EnzoInitialTurbulence
      (cycle,time,
       enzo_config->initial_turbulence_density,
       enzo_config->initial_turbulence_pressure,
       enzo_config->initial_turbulence_temperature,
       enzo::fluid_props()->gamma());
  } else if (type == "pm") {
    std::string param_str = "Initial:" + config->initial_list[index] + ":mask";
    initial = new EnzoInitialPm
      (parameters, param_str,
       cycle,time,
       enzo_config->initial_pm_field,
       enzo_config->initial_pm_mpp,
       enzo_config->initial_pm_level);
  } else if (type == "ppml_test") {
    initial = new EnzoInitialPpmlTest (cycle,time,enzo_config);
  } else if (type == "shock_tube") {
    initial = new EnzoInitialShockTube(cycle, time, p_group);
  } else if (type == "soup") {
    const int rank = enzo_config->initial_soup_rank;
    initial = new EnzoInitialSoup
      (cycle, time,
       enzo_config->initial_soup_file,
       rank,
       enzo_config->initial_soup_rotate,
       enzo_config->initial_soup_array[0],
       enzo_config->initial_soup_array[1],
       enzo_config->initial_soup_array[2],
       enzo_config->initial_soup_d_pos[0],
       enzo_config->initial_soup_d_pos[1],
       enzo_config->initial_soup_d_pos[2],
       enzo_config->initial_soup_d_size[0],
       enzo_config->initial_soup_d_size[1],
       enzo_config->initial_soup_d_size[2],
       enzo_config->initial_soup_density,
       enzo_config->initial_soup_pressure_in,
       enzo_config->initial_soup_pressure_out);
  } else if (type == "burkertbodenheimer") {
    initial = new EnzoInitialBurkertBodenheimer
      (cycle,time,
       enzo_config->initial_burkertbodenheimer_rank,
       enzo_config->initial_burkertbodenheimer_array,
       enzo_config->initial_burkertbodenheimer_radius_relative,
       enzo_config->initial_burkertbodenheimer_particle_ratio,
       enzo_config->initial_burkertbodenheimer_mass,
       enzo_config->initial_burkertbodenheimer_temperature,
       enzo_config->initial_burkertbodenheimer_densityprofile);
  } else if (type == "isolated_galaxy") {
    initial = new EnzoInitialIsolatedGalaxy (enzo_config);
  } else if (type == "merge_sinks_test") {
    initial = new EnzoInitialMergeSinksTest (enzo_config);
  } else if (type == "accretion_test") {
    initial = new EnzoInitialAccretionTest
      (cycle, time,
       enzo_config->initial_accretion_test_sink_position,
       enzo_config->initial_accretion_test_sink_velocity,
       enzo_config->initial_accretion_test_sink_mass,
       enzo_config->initial_accretion_test_gas_density,
       enzo_config->initial_accretion_test_gas_pressure,
       enzo_config->initial_accretion_test_gas_radial_velocity);
  } else if (type == "shu_collapse") {
    initial = new EnzoInitialShuCollapse
      (cycle, time,
       enzo_config->initial_shu_collapse_center,
       enzo_config->initial_shu_collapse_drift_velocity,
       enzo_config->initial_shu_collapse_truncation_radius,
       enzo_config->initial_shu_collapse_nominal_sound_speed,
       enzo_config->initial_shu_collapse_instability_parameter,
       enzo_config->initial_shu_collapse_external_density,
       enzo_config->initial_shu_collapse_central_sink_exists,
       enzo_config->initial_shu_collapse_central_sink_mass);
  } else if (type == "bb_test") {
    initial = new EnzoInitialBBTest
      (cycle, time,
       enzo_config->initial_bb_test_center,
       enzo_config->initial_bb_test_drift_velocity,
       enzo_config->initial_bb_test_mean_density,
       enzo_config->initial_bb_test_fluctuation_amplitude,
       enzo_config->initial_bb_test_truncation_radius,
       enzo_config->initial_bb_test_nominal_sound_speed,
       enzo_config->initial_bb_test_angular_rotation_velocity,
       enzo_config->initial_bb_test_external_density);
  } else {
    initial = Problem::create_initial_
      (type,index,config,parameters);
  }

  return initial;

}

//----------------------------------------------------------------------

Stopping * EnzoProblem::create_stopping_
( std::string  type, Config * config ) throw ()
/// @param type   Type of the stopping criterion to create (ignored)
/// @param config  Configuration parameter class
{
  const EnzoConfig * enzo_config = enzo::config();
  return new EnzoStopping(enzo_config->stopping_cycle,
			  enzo_config->stopping_time,
			  enzo_config->stopping_seconds,
			  enzo_config->stopping_redshift);
}

//----------------------------------------------------------------------

Refine * EnzoProblem::create_refine_
(
 std::string        type,
 int                index,
 Config *           config,
 Parameters *       parameters
 ) throw ()
{

  const EnzoConfig * enzo_config = enzo::config();

  if (type == "shock") {

    return new EnzoRefineShock
      (config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_min_refine2[index],
       config->adapt_max_coarsen2[index],
       enzo::fluid_props()->gamma(),
       enzo_config->physics_cosmology,
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index]);

  } else if (type == "particle_mass") {

    return new EnzoRefineParticleMass
      (config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index],
       config->adapt_level_exponent[index] );

  } else if (type == "mass") {

    return new EnzoRefineMass
      (config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index],
       config->adapt_field_list[index][0],
       enzo_config->adapt_mass_type[index],
       config->adapt_level_exponent[index] );

  } else {
    return Problem::create_refine_(type,index,config,parameters);
  }
}

//----------------------------------------------------------------------

Solver * EnzoProblem::create_solver_
( std::string  solver_type,
  int index_solver,
  Config * config) throw ()
/// @param solver_type   Name of the solver to create
/// @param config Configuration parameters class
{
  const EnzoConfig * enzo_config = enzo::config();

  Solver * solver = NULL;

  // Set solve type if not default "on_leaves" (solve_leaf)

  std::string solve_type_name=enzo_config->solver_solve_type[index_solver];

  int solve_type;

  if (solve_type_name=="leaf") {
    solve_type = solve_leaf;
  } else if (solve_type_name=="level") {
    solve_type = solve_level;
  } else if (solve_type_name=="tree") {
    solve_type = solve_tree;
  } else if (solve_type_name=="block") {
    solve_type = solve_block;
  } else {
    solve_type = solve_unknown;
  }

  Prolong * prolong = create_prolong_
    (enzo_config->solver_prolong[index_solver],config);
  Restrict * restrict = create_restrict_
    (enzo_config->solver_restrict[index_solver],config);

  const int index_prolong = prolong_list_.size();
  const int index_restrict = restrict_list_.size();
  prolong_list_.push_back(prolong);
  restrict_list_.push_back(restrict);

  if (solver_type == "cg") {

    solver = new EnzoSolverCg
      (enzo_config->solver_list[index_solver],
       enzo_config->solver_field_x[index_solver],
       enzo_config->solver_field_b[index_solver],
       enzo_config->solver_monitor_iter[index_solver],
       enzo_config->solver_restart_cycle[index_solver],
       solve_type,
       index_prolong,
       index_restrict,
       enzo_config->solver_min_level[index_solver],
       enzo_config->solver_max_level[index_solver],
       enzo_config->solver_iter_max[index_solver],
       enzo_config->solver_res_tol[index_solver],
       enzo_config->solver_precondition[index_solver]);

  } else if (solver_type == "dd") {

    solver = new EnzoSolverDd
      (enzo_config->solver_list[index_solver],
       enzo_config->solver_field_x[index_solver],
       enzo_config->solver_field_b[index_solver],
       enzo_config->solver_monitor_iter[index_solver],
       enzo_config->solver_restart_cycle[index_solver],
       solve_type,
       index_prolong,
       index_restrict,
       enzo_config->solver_min_level[index_solver],
       enzo_config->solver_max_level[index_solver],
       enzo_config->solver_coarse_solve[index_solver],
       enzo_config->solver_domain_solve[index_solver],
       enzo_config->solver_last_smooth[index_solver],
       enzo_config->solver_coarse_level[index_solver]);

  } else if (solver_type == "bicgstab") {

    solver = new EnzoSolverBiCgStab
      (enzo_config->solver_list[index_solver],
       enzo_config->solver_field_x[index_solver],
       enzo_config->solver_field_b[index_solver],
       enzo_config->solver_monitor_iter[index_solver],
       enzo_config->solver_restart_cycle[index_solver],
       solve_type,
       index_prolong,
       index_restrict,
       enzo_config->solver_min_level[index_solver],
       enzo_config->solver_max_level[index_solver],
       enzo_config->solver_iter_max[index_solver],
       enzo_config->solver_res_tol[index_solver],
       enzo_config->solver_precondition[index_solver],
       enzo_config->solver_coarse_level[index_solver]);

  } else if (solver_type == "diagonal") {

    solver = new EnzoSolverDiagonal
      (enzo_config->solver_list[index_solver],
       enzo_config->solver_field_x[index_solver],
       enzo_config->solver_field_b[index_solver],
       enzo_config->solver_monitor_iter[index_solver],
       enzo_config->solver_restart_cycle[index_solver],
       solve_type,
       index_prolong,
       index_restrict);

  } else if (solver_type == "jacobi") {

    solver = new EnzoSolverJacobi
      (enzo_config->solver_list[index_solver],
       enzo_config->solver_field_x[index_solver],
       enzo_config->solver_field_b[index_solver],
       enzo_config->solver_monitor_iter[index_solver],
       enzo_config->solver_restart_cycle[index_solver],
       solve_type,
       index_prolong,
       index_restrict,
       enzo_config->solver_weight[index_solver],
       enzo_config->solver_iter_max[index_solver]);

  } else if (solver_type == "mg0") {

    solver = new EnzoSolverMg0
      (enzo_config->solver_list[index_solver],
       enzo_config->solver_field_x[index_solver],
       enzo_config->solver_field_b[index_solver],
       enzo_config->solver_monitor_iter[index_solver],
       enzo_config->solver_restart_cycle[index_solver],
       solve_type,
       index_prolong,
       index_restrict,
       enzo_config->solver_min_level[index_solver],
       enzo_config->solver_max_level[index_solver],
       enzo_config->solver_iter_max[index_solver],
       enzo_config->solver_res_tol[index_solver],
       enzo_config->solver_pre_smooth[index_solver],
       enzo_config->solver_coarse_solve[index_solver],
       enzo_config->solver_post_smooth[index_solver],
       enzo_config->solver_last_smooth[index_solver],
       enzo_config->solver_coarse_level[index_solver]);

  } else {
    // Not an Enzo Solver--try base class Cello Solver
    solver = Problem::create_solver_ (solver_type, index_solver,config);
  }

  ASSERT1 ("EnzoProblem::create_solver()",
	   "Unknown solver %s",
	   solver_type.c_str(),
	   solver != NULL);

  solver->set_index(index_solver);

  return solver;
}

//----------------------------------------------------------------------

Compute * EnzoProblem::create_compute
( std::string name,
  Config * config ) throw()
/// @param name  Name of the compute to create
{

  Compute * compute = 0;

  TRACE1("EnzoProblem::create_compute %s",name.c_str());

  const EnzoConfig * enzo_config = enzo::config();

  if (name == "temperature") {

    compute = new EnzoComputeTemperature(enzo::fluid_props(),
                                         enzo_config->physics_cosmology);

  } else if (name == "pressure"){

    compute = new EnzoComputePressure(enzo::fluid_props()->gamma(),
                                      enzo_config->physics_cosmology);

#ifdef CONFIG_USE_GRACKLE
  } else if (name == "cooling_time"){

    compute = new EnzoComputeCoolingTime();

#endif
  } else {

    // Fallback to Cello method's
    compute = Problem::create_compute (name,config);

    ASSERT2("EnzoProblem::create_compute",
            "Compute created %s does not match compute requested %s",
            compute->name().c_str(),name.c_str(),
            compute->name() == name);
  }

  return compute;
}


//----------------------------------------------------------------------

Method * EnzoProblem::create_method_
( std::string  name,
  int index_method,
  Config * config,
  const Factory * factory) throw ()
/// @param name   Name of the method to create
/// @param config Configuration parameters class
{
  Method * method = 0;

  // historically, this method would always call method->set_courant after
  // building a new method object. But, with this new p_group approach, each
  // method objects should opt out of this approach (so that they can set
  // appropriate default courant factors)
  // - in cases where the courant factor is not used, we may not explicitly set
  //   this variable to true (since it doesn't really matter)
  bool skip_auto_courant = false;

  // move creation of p_group up the call stack?
  ASSERT("Problem::create_method_", "Something is wrong", cello::simulation());
  Parameters* parameters = cello::simulation()->parameters();
  const std::string root_path =
    ("Method:" + parameters->list_value_string(index_method, "Method:list"));
  ParameterGroup p_group(*parameters, root_path);

  const EnzoConfig * enzo_config = enzo::config();

  // The following 2 lines may need to be updated in the future
  const std::vector<std::string>& mlist = enzo_config->method_list;
  const bool store_fluxes_for_corrections =
    std::find(mlist.begin(), mlist.end(), "flux_correct") != mlist.end();

  TRACE1("EnzoProblem::create_method %s",name.c_str());
  if (name == "ppm") {
    method = new EnzoMethodPpm
      (store_fluxes_for_corrections,
       enzo_config->ppm_diffusion,
       enzo_config->ppm_flattening,
       enzo_config->ppm_pressure_free,
       enzo_config->ppm_steepening,
       enzo_config->ppm_use_minimum_pressure_support,
       enzo_config->ppm_minimum_pressure_support_parameter);

  } else if (name == "ppml") {

    method = new EnzoMethodPpml(p_group);
    skip_auto_courant = true;

  } else if (name == "pm_deposit") {

    method = new EnzoMethodPmDeposit (p_group);

  } else if (name == "pm_update") {

    method = new EnzoMethodPmUpdate(p_group);

  } else if (name == "heat") {

    method = new EnzoMethodHeat(p_group);
    skip_auto_courant = true;

#ifdef CONFIG_USE_GRACKLE

  } else if (name == "grackle") {

    method = new EnzoMethodGrackle
      (enzo_config->method_grackle_chemistry,
       enzo_config->method_grackle_use_cooling_timestep,
       enzo_config->method_grackle_radiation_redshift,
       enzo_config->physics_cosmology_initial_redshift,
       enzo::simulation()->time());

#endif /* CONFIG_USE_GRACKLE */

  } else if (name == "balance") {

    method = new EnzoMethodBalance;

  } else if (name == "turbulence") {

    method = new EnzoMethodTurbulence
      (enzo_config->method_turbulence_edot,
       enzo_config->initial_turbulence_density,
       enzo_config->initial_turbulence_temperature,
       enzo_config->method_turbulence_mach_number,
       enzo_config->physics_cosmology);

  } else if (name == "cosmology") {

    method = new EnzoMethodCosmology;

  } else if (name == "comoving_expansion") {

    bool comoving_coordinates = enzo_config->physics_cosmology;

    method = new EnzoMethodComovingExpansion ( comoving_coordinates );

  } else if (name == "gravity") {

    std::string solver_name = enzo_config->method_gravity_solver;

    int index_solver = enzo_config->solver_index.at(solver_name);

    ASSERT1 ("EnzoProblem::create_solver_()",
	     "Cannot find solver \"%s\"",
	     solver_name.c_str(),
	     0 <= index_solver && index_solver < enzo_config->num_solvers);

    Prolong * prolong = create_prolong_
      (config->method_prolong[index_method],config);

    const int index_prolong = prolong_list_.size();
    prolong_list_.push_back(prolong);

    method = new EnzoMethodGravity
      (
       enzo_config->solver_index.at(solver_name),
       enzo_config->method_gravity_order,
       enzo_config->method_gravity_accumulate,
       index_prolong,
       enzo_config->method_gravity_dt_max);

  } else if (name == "mhd_vlct") {

    method = new EnzoMethodMHDVlct
      (enzo_config->method_vlct_riemann_solver,
       enzo_config->method_vlct_time_scheme,
       enzo_config->method_vlct_reconstruct_method,
       enzo_config->method_vlct_theta_limiter,
       enzo_config->method_vlct_mhd_choice,
       store_fluxes_for_corrections);

  } else if (name == "background_acceleration") {

    // If self-gravity is calculated, we do not need to zero
    // out the acceleration field from the previous time step
    // before adding the background accelerations
    bool zero_acceleration = true;
    for (int index = 0; index < method_list_.size(); index++){
      if (method_list_[index]->name() == "gravity"){
        zero_acceleration = false;
        break;
      }
    }

    method = new EnzoMethodBackgroundAcceleration
      (zero_acceleration);

  } else if (name == "star_maker") {

    // should generalize this to enable multiple maker types
    if (enzo_config->method_star_maker_flavor == "stochastic"){
      method = new EnzoMethodStarMakerStochasticSF();
    } else if (enzo_config->method_star_maker_flavor == "STARSS" ||
               enzo_config->method_star_maker_flavor == "starss") {
      method = new EnzoMethodStarMakerSTARSS();
    } else{ // does not do anything
      method = new EnzoMethodStarMaker();
    }

  } else if (name == "feedback") {

    // need a similar type swtich as in star maker
    if (enzo_config->method_feedback_flavor == "distributed"){
      method = new EnzoMethodDistributedFeedback();
    } else if (enzo_config->method_feedback_flavor == "STARSS" ||
               enzo_config->method_feedback_flavor == "starss") {
      method = new EnzoMethodFeedbackSTARSS();
    }  else { // does not do anything
      method = new EnzoMethodFeedback();
    }

  } else if (name == "m1_closure") {

    method = new EnzoMethodM1Closure(enzo_config->method_m1_closure_N_groups);

  } else if (name == "check") {

    // Method for checkpointing the simulation
    method = new EnzoMethodCheck
      (enzo_config->method_check_num_files,
       enzo_config->method_check_ordering,
       enzo_config->method_check_dir,
       enzo_config->method_check_monitor_iter,
       enzo_config->method_check_include_ghosts);

  } else if (name == "merge_sinks") {

    method = new EnzoMethodMergeSinks
      (enzo_config->method_merge_sinks_merging_radius_cells);

  } else if (name == "accretion") {

    if (enzo_config->method_accretion_flavor == "threshold") {
      method = new EnzoMethodThresholdAccretion
        (
         enzo_config->method_accretion_accretion_radius_cells,
         enzo_config->method_accretion_physical_density_threshold_cgs,
         enzo_config->method_accretion_max_mass_fraction
         );
    } else if (enzo_config->method_accretion_flavor == "bondi_hoyle") {
      method = new EnzoMethodBondiHoyleAccretion
        (
         enzo_config->method_accretion_accretion_radius_cells,
         enzo_config->method_accretion_physical_density_threshold_cgs,
         enzo_config->method_accretion_max_mass_fraction
         );
    } else if (enzo_config->method_accretion_flavor == "flux") {
      method = new EnzoMethodFluxAccretion
        (
         enzo_config->method_accretion_accretion_radius_cells,
         enzo_config->method_accretion_physical_density_threshold_cgs,
         enzo_config->method_accretion_max_mass_fraction
         );
    } else if (enzo_config->method_accretion_flavor == "dummy"){
      method = new EnzoMethodAccretion
        (
         enzo_config->method_accretion_accretion_radius_cells,
         enzo_config->method_accretion_physical_density_threshold_cgs,
         enzo_config->method_accretion_max_mass_fraction
         );
    } else {
      ERROR1("EnzoProblem::create_method_",
             "\"accretion\" method has flavor \"%s\", which is not one of the possible options: "
	     "\"threshold\", \"bondi_hoyle\", \"flux\", or \"dummy\"",
	     enzo_config->method_accretion_flavor.c_str());
    }
  } else if (name == "sink_maker") {

    method = new EnzoMethodSinkMaker(p_group);

  } else {

    // Fallback to Cello method's
    method = Problem::create_method_ (name, index_method,config,
                                      enzo::simulation()->factory());

  }

  if (method) {

    // set the method's courant safety factor
    if (!skip_auto_courant){
      method->set_courant(config->method_courant[index_method]);
    }

    ASSERT2("EnzoProblem::create_method",
	    "Method created %s does not match method requested %s",
	    method->name().c_str(),name.c_str(),
	    method->name() == name);
  }

  return method;
}

//----------------------------------------------------------------------

Prolong * EnzoProblem::create_prolong_
( std::string  type,
  Config *     config ) throw ()
{

  Prolong * prolong = 0;

  const EnzoConfig * enzo_config = enzo::config();

  if (type == "enzo") {
    prolong = new EnzoProlong
      (enzo_config->prolong_enzo_type,
       enzo_config->prolong_enzo_positive,
       enzo_config->prolong_enzo_use_linear);
  } else {
    prolong = Problem::create_prolong_(type,config);
  }

  return prolong;

}

//----------------------------------------------------------------------

Physics * EnzoProblem::create_physics_
( std::string  type,
   int index,
   Config * config,
   Parameters * parameters) throw ()
{

  Physics * physics = NULL;
  const EnzoConfig * enzo_config = enzo::config();

  if (type == "cosmology") {

    physics = new EnzoPhysicsCosmology
      (
       enzo_config->physics_cosmology_hubble_constant_now,
       enzo_config->physics_cosmology_omega_matter_now,
       enzo_config->physics_cosmology_omega_baryon_now,
       enzo_config->physics_cosmology_omega_cdm_now,
       enzo_config->physics_cosmology_omega_lamda_now,
       enzo_config->physics_cosmology_comoving_box_size,
       enzo_config->physics_cosmology_max_expansion_rate,
       enzo_config->physics_cosmology_initial_redshift,
       enzo_config->physics_cosmology_final_redshift
       );

  } else if (type == "fluid_props") {

    physics = new EnzoPhysicsFluidProps
      (
       enzo_config->physics_fluid_props_de_config,
       enzo_config->physics_fluid_props_fluid_floor_config,
       enzo_config->physics_fluid_props_eos_variant,
       enzo_config->physics_fluid_props_mol_weight
       );

  } else if (type == "gravity") {

    // note: it may make sense to convert EnzoProblem::initialize_physics into
    // a virtual method and provide a custom implementation that reorders the
    // physics object initialization to avoid the following problem
    // - unlike things such as methods, initializers, boundaries, or refinement
    //   criteria, ordering of physics object shouldn't matter to a user
    // - if we did that, we could consolidate that method with the
    //   initialize_physics_coda_ method

    for (std::size_t i = index; i < enzo_config->physics_list.size(); i++) {
      ASSERT("EnzoProblem::create_physics_",
             "a \"cosmology\" physics object MUST NOT follow a \"gravity\" "
             "object (it's okay if it comes before the \"gravity\" object)",
             enzo_config->physics_list[i] != "cosmology");
    }
    physics = new EnzoPhysicsGravity
      (enzo_config->physics_gravity_grav_constant_codeU);

  } else {

    physics = Problem::create_physics_
      (type,index,config,parameters);
  }

  return physics;

}

//----------------------------------------------------------------------

void EnzoProblem::initialize_physics_coda_(Config * config,
                                           Parameters * parameters) throw()
{
  // if EnzoPhysicsFluidProps or EnzoPhysicsGravity don't already exist,
  // initialize them (this is required for backwards compatability)
  const std::vector<std::string> required = {"fluid_props", "gravity"};
  for (const std::string& name: required) {
    if (physics(name) != nullptr) { continue; }
    physics_list_.push_back(create_physics_(name, physics_list_.size(),
                                            config, parameters));
  }

  // in the future, we might want to move the following snippet from
  // EnzoInitialCosmology::EnzoInitialCosmology to this function:
  //  EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology *)
  //  this->physics("cosmology");
  //  if (cosmology != nullptr) { enzo::units()->set_cosmology(cosmology); }
  // Doing this could resolve some issues encountered in EnzoMethodGrackle more
  // elegantly that the existing work-around
}

//----------------------------------------------------------------------

Units * EnzoProblem::create_units_ (  Config * config  ) throw ()
{
  EnzoUnits * units = new EnzoUnits;

  if (config->units_mass == 1.0) {

    units->set_using_density (config->units_length,
			      config->units_density,
			      config->units_time);

  } else if (config->units_density == 1.0) {

    units->set_using_mass (config->units_length,
			   config->units_mass,
			   config->units_time);
  } else {

    ERROR("Problem::create_units_",
	  "Cannot set both Units:density and Units:mass parameters");
  }

  return units;

}

//----------------------------------------------------------------------

Restrict * EnzoProblem::create_restrict_
(
 std::string  type,
 Config * config ) throw ()
{

  Restrict * restrict = 0;

  restrict = Problem::create_restrict_(type,config);

  return restrict;

}

//----------------------------------------------------------------------
