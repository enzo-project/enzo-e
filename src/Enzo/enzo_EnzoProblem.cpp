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
  } else if (type == "inflow") {
    boundary = new EnzoBoundary (axis,face,mask,boundary_type_inflow);
  } else if (type == "periodic") {
    boundary = new EnzoBoundary (axis,face,mask,boundary_type_periodic);
  } else {
    boundary = Problem::create_boundary_(type,index,config,parameters);
  }

  return boundary;
}

// //----------------------------------------------------------------------

// Timestep * EnzoProblem::create_timestep_
// (
//  std::string  type,
//  Config * config
//  ) throw ()
// /// @param type   Type of the timestep method to create
// {
//   // Timestep * timestep = 0;

//   // if (type == "ppml") {
//   //   timestep = new EnzoTimestepPpml;
//   // } else if (type == "ppm") {
//   //   timestep = new EnzoTimestep;
//   // } else {
//   //   timestep = Problem::create_timestep_(type,config);
//   // }

//   return NULL;
// }

// //----------------------------------------------------------------------

Initial * EnzoProblem::create_initial_ 
(
 std::string  type,
 Parameters * parameters,
 Config * config,
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

  if (type == "implosion_2d") {
    initial = new EnzoInitialImplosion2(cycle,time);
  } else if (type == "sedov_array_2d") {
    initial = new EnzoInitialSedovArray2(static_cast<EnzoConfig *>(config));
  } else if (type == "sedov_array_3d") {
    initial = new EnzoInitialSedovArray3(static_cast<EnzoConfig *>(config));
  } else {
    initial = Problem::create_initial_
      (type,config,parameters,field_descr,group_process);
  }

  return initial;

}

//----------------------------------------------------------------------

Method * EnzoProblem::create_method_ 
( std::string  type,  Config * config) throw ()
/// @param type   Type of the method to create
{

  Method * method = 0;

  EnzoConfig * enzo_config = (EnzoConfig *) config;

  TRACE1("EnzoProblem::create_method %s",type.c_str());
  if (type == "ppm") {
    method = new EnzoMethodPpm;
  } else if (type == "ppml") {
    method = new EnzoMethodPpml;
  } else if (type == "heat") {
    method = new EnzoMethodHeat(enzo_config->enzo_method_heat_alpha,
				enzo_config->field_courant);
  } else {
    method = Problem::create_method_(type,config);
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
      (static_cast<EnzoConfig *>(config)->enzo_interpolation_method);

  } else if (type == "MC1") {
    
    prolong = new EnzoProlongMC1
      (static_cast<EnzoConfig *>(config)->enzo_interpolation_method);

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
      (static_cast<EnzoConfig *>(config)->enzo_interpolation_method);

  } else {

    restrict = Problem::create_restrict_(type,config);

  }

  return restrict;
  
}

//----------------------------------------------------------------------

