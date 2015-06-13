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

  Mask * mask = 0;
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

  if (type == "implosion_2d") {
    initial = new EnzoInitialImplosion2(cycle,time);
  } else if (type == "sedov_array_2d") {
    initial = new EnzoInitialSedovArray2(enzo_config);
  } else if (type == "sedov_array_3d") {
    initial = new EnzoInitialSedovArray3(enzo_config);
#ifdef CONFIG_USE_GRACKLE
  } else if (type == "grackle_test") {
    initial = new EnzoInitialGrackleTest(enzo_config);
#endif /* CONFIG_USE_GRACKLE */
  } else if (type == "turbulence") {
    initial = new EnzoInitialTurbulence 
      (cycle,time, 
       enzo_config->initial_turbulence_density,
       enzo_config->initial_turbulence_pressure,
       enzo_config->initial_turbulence_temperature,
       enzo_config->field_gamma);
  } else {
    initial = Problem::create_initial_
      (type,config,parameters,field_descr);
  }

  return initial;

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
       config->mesh_min_refine[index],
       config->mesh_max_coarsen[index],
       config->mesh_min_refine2[index],
       config->mesh_max_coarsen2[index],
       enzo_config->field_gamma,
       enzo_config->physics_cosmology,
       config->mesh_refine_output[index]);
  } else {
    return Problem::create_refine_(type,config,parameters,field_descr,index);
  }
}


//----------------------------------------------------------------------

Method * EnzoProblem::create_method_ 
( std::string  name,  
  Config * config,
  const FieldDescr * field_descr) throw ()
/// @param name   Name of the method to create
/// @param config Configuration parameters class
/// @param field_descr Field descriptor
{

  Method * method = 0;

  EnzoConfig * enzo_config = static_cast<EnzoConfig *>(config);

  TRACE1("EnzoProblem::create_method %s",name.c_str());
  if (name == "ppm") {
    method = new EnzoMethodPpm (enzo_config);
  } else if (name == "ppml") {
    method = new EnzoMethodPpml(enzo_config);
  } else if (name == "heat") {
    method = new EnzoMethodHeat
      (field_descr,
       enzo_config->method_heat_alpha,
       enzo_config->field_courant);
  } else if (name == "null") {
    method = new EnzoMethodNull
      (enzo_config->method_null_dt);
#ifdef CONFIG_USE_GRACKLE
  } else if (name == "grackle") {
    method = new EnzoMethodGrackle (enzo_config);
#endif /* CONFIG_USE_GRACKLE */
  } else if (name == "turbulence") {
    method = new EnzoMethodTurbulence 
      (enzo_config->method_turbulence_edot,
       enzo_config->initial_turbulence_density,
       enzo_config->initial_turbulence_temperature,
       enzo_config->method_turbulence_mach_number,
       enzo_config->physics_cosmology);
  } else if (name == "gravity_cg") {
    const bool is_singular = is_periodic();
    int rank = config->mesh_root_rank;
    method = new EnzoMethodGravityCg
      (field_descr, rank,
       enzo_config->method_gravity_cg_grav_const,
       enzo_config->method_gravity_cg_iter_max,
       enzo_config->method_gravity_cg_res_tol,
       enzo_config->method_gravity_cg_monitor_iter,
       is_singular,
       enzo_config->method_gravity_cg_diag_precon );
  } else if (name == "gravity_mg") {
    const bool is_singular = is_periodic();
    int rank = config->mesh_root_rank;
    Restrict * restrict = 
      create_restrict_(enzo_config->method_gravity_mg_restrict,config);
    Prolong * prolong = 
      create_prolong_(enzo_config->method_gravity_mg_prolong,config);
    Compute * smooth = NULL;
    if (enzo_config->method_gravity_mg_smooth == "jacobi") {
      smooth = new EnzoComputeSmoothJacobi 
	("X","R","D",
	 enzo_config->method_gravity_mg_smooth_weight,
	 field_descr);
    }
    std::string type = enzo_config->method_gravity_mg_type;
    if (type == "mlat") {
      method = new EnzoMethodGravityMlat
	(field_descr, rank,
	 enzo_config->method_gravity_mg_grav_const,
	 enzo_config->method_gravity_mg_iter_max,
	 enzo_config->method_gravity_mg_res_tol,
	 enzo_config->method_gravity_mg_monitor_iter,
	 is_singular,  smooth, restrict,  prolong,
	 enzo_config->method_gravity_mg_level_min,
	 enzo_config->method_gravity_mg_level_max);
    } else if (type == "mg0") {
      method = new EnzoMethodGravityMg0
	(field_descr, rank,
	 enzo_config->method_gravity_mg_grav_const,
	 enzo_config->method_gravity_mg_iter_max,
	 enzo_config->method_gravity_mg_monitor_iter,
	 is_singular,  smooth, restrict,  prolong,
	 enzo_config->method_gravity_mg_level_min,
	 enzo_config->method_gravity_mg_level_max);
    } else {
      ERROR1 ("EnzoProblem::create_method",
	       "Unknown gravity_mg type %s",
	       type.c_str());
	       
    }
  } else if (name == "gravity_bicgstab") {
    method = new EnzoMethodGravityBiCGStab
      (field_descr,
       enzo_config->method_gravity_bicgstab_iter_max,
       enzo_config->method_gravity_bicgstab_res_tol);
  } else {
    method = Problem::create_method_ (name,config, field_descr);
  }

  ASSERT2("EnzoProblem::create_method",
	  "Method created %s does not match method requested %s",
	  method->name().c_str(),name.c_str(),
	  method->name() == name);

  return method;
}

//----------------------------------------------------------------------

Prolong * EnzoProblem::create_prolong_ 
( std::string  type,
  Config *     config ) throw ()
{

  Prolong * prolong = 0;

  EnzoConfig * enzo_config = static_cast<EnzoConfig *>(config);

  if (type == "enzo") {
    
    prolong = new EnzoProlong (enzo_config->interpolation_method);

  } else if (type == "MC1") {
    
    prolong = new EnzoProlongMC1 (enzo_config->interpolation_method);

  } else if (type == "poisson") {
    
    prolong = new EnzoProlongPoisson;

  } else {

    prolong = Problem::create_prolong_(type,config);

  }

  return prolong;
  
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

