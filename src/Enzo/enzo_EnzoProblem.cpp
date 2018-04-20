// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProblem.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-03
/// @brief    Implementation of EnzoProblem class
///
/// 

#include "enzo.hpp"

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
 Parameters * parameters,
 const FieldDescr * field_descr
 ) throw ()
{
  
  //--------------------------------------------------
  // parameter: Initial : cycle
  // parameter: Initial : time
  //--------------------------------------------------

  Initial * initial = 0;

  int cycle   = config->initial_cycle;
  double time = config->initial_time;

  EnzoConfig * enzo_config = static_cast<EnzoConfig *>(config);

  if (type == "music") {
    
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
       enzo_config->field_gamma,
       enzo_config->initial_cosmology_temperature
       );
  } else if (type == "turbulence") {
    initial = new EnzoInitialTurbulence 
      (cycle,time, 
       enzo_config->initial_turbulence_density,
       enzo_config->initial_turbulence_pressure,
       enzo_config->initial_turbulence_temperature,
       enzo_config->field_gamma);
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
  } else {
    initial = Problem::create_initial_
      (type,index,config,parameters,field_descr);
  }

  return initial;

}

//----------------------------------------------------------------------

Stopping * EnzoProblem::create_stopping_ 
( std::string  type, Config * config ) throw ()
/// @param type   Type of the stopping criterion to create (ignored)
/// @param config  Configuration parameter class
{
  EnzoConfig * enzo_config = static_cast<EnzoConfig *>(config);
  return new EnzoStopping(enzo_config->stopping_cycle,
			  enzo_config->stopping_time,
			  enzo_config->stopping_seconds,
			  enzo_config->stopping_redshift);
}

//----------------------------------------------------------------------

Refine * EnzoProblem::create_refine_
(
 std::string        type,
 Config *           config,
 Parameters *       parameters,
 const FieldDescr * field_descr,
 int                index
 ) throw ()
{ 

  EnzoConfig * enzo_config = static_cast<EnzoConfig *>(config);

  if (type == "shock") {

    return new EnzoRefineShock 
      (field_descr,
       config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_min_refine2[index],
       config->adapt_max_coarsen2[index],
       enzo_config->field_gamma,
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
    return Problem::create_refine_(type,config,parameters,field_descr,index);
  }
}

//----------------------------------------------------------------------

Solver * EnzoProblem::create_solver_ 
( std::string  solver_type,  
  Config * config,
  int index_solver,
  FieldDescr * field_descr,
  const ParticleDescr * particle_descr) throw ()
/// @param solver_type   Name of the solver to create
/// @param config Configuration parameters class
/// @param field_descr Field descriptor
{
  EnzoConfig * enzo_config = static_cast<EnzoConfig *>(config);

  Solver * solver = NULL;
  
  int rank = config->mesh_root_rank;

  if (solver_type == "cg") {

    solver = new EnzoSolverCg
      (field_descr,
       enzo_config->solver_monitor_iter [index_solver],
       enzo_config->solver_restart_cycle [index_solver],
       rank,
       enzo_config->solver_iter_max     [index_solver],
       enzo_config->solver_res_tol      [index_solver],
       enzo_config->solver_min_level    [index_solver],
       enzo_config->solver_max_level    [index_solver],
       enzo_config->solver_precondition [index_solver],
       enzo_config->solver_local        [index_solver]
       );

  } else if (solver_type == "bicgstab") {

    solver = new EnzoSolverBiCgStab
      (field_descr,
       enzo_config->solver_monitor_iter[index_solver],
       enzo_config->solver_restart_cycle [index_solver],
       rank,
       enzo_config->solver_iter_max[index_solver],
       enzo_config->solver_res_tol[index_solver],
       enzo_config->solver_min_level[index_solver],
       enzo_config->solver_max_level[index_solver],
       enzo_config->solver_precondition[index_solver]) ;

  } else if (solver_type == "diagonal") {

    solver = new EnzoSolverDiagonal;

  } else if (solver_type == "jacobi") {

    solver = new EnzoSolverJacobi
      (field_descr,
       enzo_config->solver_weight[index_solver],
       enzo_config->solver_iter_max[index_solver]);

  } else if (solver_type == "mg0") {

    Restrict * restrict = 
      create_restrict_(enzo_config->solver_restrict[index_solver],config);
    Prolong * prolong = 
      create_prolong_(enzo_config->solver_prolong[index_solver],config);

    solver = new EnzoSolverMg0
      (field_descr,
       enzo_config->solver_monitor_iter[index_solver],
       enzo_config->solver_restart_cycle [index_solver],
       rank,
       enzo_config->solver_iter_max[index_solver],
       enzo_config->solver_res_tol[index_solver],
       enzo_config->solver_pre_smooth[index_solver],
       enzo_config->solver_coarse_solve[index_solver],
       enzo_config->solver_post_smooth[index_solver],
       enzo_config->solver_last_smooth[index_solver],
       restrict,  prolong,
       enzo_config->solver_min_level[index_solver],
       enzo_config->solver_max_level[index_solver]);

  } else {
    // Not an Enzo Solver--try base class Cello Solver
    solver = Problem::create_solver_ 
      (solver_type,config, index_solver,(FieldDescr *)field_descr,particle_descr);
    
  }

  ASSERT1 ("EnzoProblem::create_solver()",
	   "Unknown solver %s",
	   solver_type.c_str(),
	   solver != NULL);

  solver->set_index(index_solver);

  return solver;
}    

//----------------------------------------------------------------------

Method * EnzoProblem::create_method_ 
( std::string  name,  
  Config * config,
  int index_method,
  const FieldDescr * field_descr,
  const ParticleDescr * particle_descr) throw ()
/// @param name   Name of the method to create
/// @param config Configuration parameters class
/// @param field_descr Field descriptor
{

  Method * method = 0;

  EnzoConfig * enzo_config = static_cast<EnzoConfig *>(config);
 
  TRACE1("EnzoProblem::create_method %s",name.c_str());
  
  if (name == "ppm") {

    method = new EnzoMethodPpm (field_descr, enzo_config);

  } else if (name == "hydro") {

    method = new EnzoMethodHydro
      (field_descr,
       enzo_config->method_hydro_method,
       enzo_config->field_gamma,
       enzo_config->physics_gravity,
       enzo_config->physics_cosmology,
       enzo_config->method_hydro_dual_energy,
       enzo_config->method_hydro_dual_energy_eta_1,
       enzo_config->method_hydro_dual_energy_eta_2,
       enzo_config->method_hydro_reconstruct_method,
       enzo_config->method_hydro_reconstruct_conservative,
       enzo_config->method_hydro_reconstruct_positive,
       enzo_config->ppm_density_floor,
       enzo_config->ppm_pressure_floor,
       enzo_config->ppm_pressure_free,
       enzo_config->ppm_diffusion,
       enzo_config->ppm_flattening,
       enzo_config->ppm_steepening,
       enzo_config->method_hydro_riemann_solver
       );

  } else if (name == "ppml") {

    method = new EnzoMethodPpml (field_descr, enzo_config);

  } else if (name == "pm_deposit") {

    method = new EnzoMethodPmDeposit  
      (field_descr, particle_descr, 
       enzo_config->method_pm_deposit_type, 0.5);

  } else if (name == "pm_update") {

    method = new EnzoMethodPmUpdate  
      (field_descr, particle_descr, enzo_config->method_pm_update_max_dt);

  } else if (name == "heat") {

    method = new EnzoMethodHeat
      (field_descr,
       enzo_config->method_heat_alpha,
       config->method_courant[index_method]);

  } else if (name == "null") {

    method = new EnzoMethodNull
      (field_descr, enzo_config->method_null_dt);

#ifdef CONFIG_USE_GRACKLE
    //--------------------------------------------------
  } else if (name == "grackle") {
    method = new EnzoMethodGrackle (enzo_config,field_descr);
#endif /* CONFIG_USE_GRACKLE */
    
  } else if (name == "turbulence") {
    
    method = new EnzoMethodTurbulence 
      (field_descr,
       enzo_config->method_turbulence_edot,
       enzo_config->initial_turbulence_density,
       enzo_config->initial_turbulence_temperature,
       enzo_config->method_turbulence_mach_number,
       enzo_config->physics_cosmology);

  } else if (name == "cosmology") {

    method = new EnzoMethodCosmology(field_descr);

  } else if (name == "comoving_expansion") {

    bool comoving_coordinates = enzo_config->physics_cosmology;
    
    method = new EnzoMethodComovingExpansion
      ( field_descr, comoving_coordinates );

  } else if (name == "gravity") {

    std::string solver_name = enzo_config->method_gravity_solver;

    int index_solver = 0;
    while (index_solver < enzo_config->num_solvers &&
	   enzo_config->solver_list[index_solver] != solver_name) {
      ++index_solver;
    }
    ASSERT1 ("EnzoProblem::create_solver_()",
	     "Cannot find solver \"%s\"",
	     solver_name.c_str(),
	     index_solver < enzo_config->num_solvers);

   // WARNING: method_gravity_accumulate is currently only implemented
   // for non-AMR problems

  method = new EnzoMethodGravity
      (field_descr, enzo_config->solver_index[solver_name],
       enzo_config->method_gravity_grav_const,
       enzo_config->method_gravity_order,
       enzo_config->method_gravity_accumulate);
      
  } else {

    // Fallback to Cello method's
    method = Problem::create_method_ 
      (name,config, index_method,field_descr,particle_descr);

  }

  if (method) {

    // set the method's courant safety factor
    method->set_courant(config->method_courant[index_method]);

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

  prolong = Problem::create_prolong_(type,config);

  return prolong;
  
}

//----------------------------------------------------------------------

Physics * EnzoProblem::create_physics_ 
( std::string  type,
   int index,
   Config * config,
   Parameters * parameters,
   const FieldDescr * field_descr) throw ()
{

  Physics * physics = NULL;

  if (type == "cosmology") {

    EnzoConfig * enzo_config = static_cast<EnzoConfig *>(config);

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

  } else {
    
    physics = Problem::create_physics_
      (type,index,config,parameters,field_descr);
    
  }

  return physics;
  
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
	  "Cannot set both Units:density and Units:time parameters");
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

  if (type == "enzo") {
    
    restrict = new EnzoRestrict 
      (static_cast<EnzoConfig *>(config)->interpolation_method);

  } else {

    restrict = Problem::create_restrict_(type,config);

  }

  return restrict;
  
}

//----------------------------------------------------------------------

