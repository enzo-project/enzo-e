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
 const FieldDescr * field_descr,
 const GroupProcess * group_process
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
      (type,config,parameters,field_descr,group_process);
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
       config->mesh_refine_output[index]);
  } else {
    return Problem::create_refine_(type,config,parameters,field_descr,index);
  }
}


//----------------------------------------------------------------------

Method * EnzoProblem::create_method_ 
( std::string  type,  
  Config * config,
  const FieldDescr * field_descr) throw ()
/// @param type   Type of the method to create
/// @param config Configuration parameters class
{

  Method * method = 0;

  EnzoConfig * enzo_config = (EnzoConfig *) config;

  TRACE1("EnzoProblem::create_method %s",type.c_str());
  if (type == "ppm") {
    method = new EnzoMethodPpm;
  } else if (type == "ppm3") {
    method = new EnzoMethodPpm3;
  } else if (type == "ppml") {
    method = new EnzoMethodPpml;
  } else if (type == "heat") {
    method = new EnzoMethodHeat
      (enzo_config->method_heat_alpha,
       enzo_config->field_courant);
  } else if (type == "null") {
    method = new EnzoMethodNull
      (enzo_config->method_null_dt);
#ifdef CONFIG_USE_GRACKLE
  } else if (type == "grackle") {
    method = new EnzoMethodGrackle (enzo_config);
#endif /* CONFIG_USE_GRACKLE */
  } else if (type == "turbulence") {
    method = new EnzoMethodTurbulence 
      (enzo_config->method_turbulence_edot,
       enzo_config->initial_turbulence_density,
       enzo_config->initial_turbulence_temperature,
       enzo_config->method_turbulence_mach_number);
  } else if (type == "pressure") {
    method = new EnzoMethodPressure (enzo_config->field_gamma);
  } else if (type == "temperature") {
    method = new EnzoMethodTemperature 
      (enzo_config->ppm_density_floor,
       enzo_config->ppm_temperature_floor,
       enzo_config->ppm_mol_weight);
  } else if (type == "gravity_bicgstab") {
    method = new EnzoMethodGravityBiCGStab
      (enzo_config->method_gravity_bicgstab_iter_max,
       enzo_config->method_gravity_bicgstab_res_tol);
  } else {
    method = Problem::create_method_ (type,config, field_descr);
  }

  return method;
}

//----------------------------------------------------------------------

Prolong * EnzoProblem::create_prolong_ 
( std::string  type,
  Config *     config ) throw ()
{

  Prolong * prolong = 0;

  if (type == "enzo") {
    
    prolong = new EnzoProlong 
      (static_cast<EnzoConfig *>(config)->interpolation_method);

  } else if (type == "MC1") {
    
    prolong = new EnzoProlongMC1
      (static_cast<EnzoConfig *>(config)->interpolation_method);

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

