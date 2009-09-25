//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */


/** 
 *********************************************************************
 *
 * @file      method_ppm.cpp
 * @brief     Implements the PPM Enzo hydrodynamics Method
 * @author    James
 * @date      Sat Jul 18 13:50:46 PDT 2009
 *
 * DESCRIPTION 
 * 
 *    PPM Enzo hydrodynamics Method
 *
 * PACKAGES
 *
 *    Field
 *    Particle
 * 
 * INCLUDES
 *  
 *    
 *
 * PUBLIC FUNCTIONS
 *  
 *    MethodPpm::MethodPpm()
 *    MethodPpm::initialize()
 *    MethodPpm::add_argument()
 *    MethodPpm::apply()
 *
 * PRIVATE FUCTIONS
 *  
 * 
 *
 * $Id$
 *
 *********************************************************************
 */

#include <string>

#include "method.hpp"
 
/**
 *********************************************************************
 *
 * @return        There is no return value
 *
 * Create a new uninitialized MethodPpm
 *
 *********************************************************************
 */

MethodPpm::MethodPpm() throw()
  : field_name_density_(),
    field_name_velocity_x_(),
    field_name_velocity_y_(),
    field_name_velocity_z_(),
    field_name_total_energy_(),
    field_name_pressure_()
{
}

/**
 *********************************************************************
 *
 * @param  method The MethodPpm being copied 
 * @return        There is no return value
 *
 * Create a new copy of the given MethodPpm
 *
 *********************************************************************
 */

MethodPpm::MethodPpm(const MethodPpm & method) throw()
{
  field_name_density_      = method.field_name_density_;
  field_name_velocity_x_   = method.field_name_velocity_x_;
  field_name_velocity_y_   = method.field_name_velocity_y_;
  field_name_velocity_z_   = method.field_name_velocity_z_;
  field_name_total_energy_ = method.field_name_total_energy_;
  field_name_pressure_     = method.field_name_pressure_;
}

/**
 *********************************************************************
 *
 * @return        There is no return value
 *
 * Initialize the MethodPpm
 *
 *********************************************************************
 */

void MethodPpm::initialize() throw()
{
}

/**
 *********************************************************************
 *
 * @param  name   Name of the Field or Particle set, or tag
 * @return        There is no return value
 *
 * Specify a Field or Particle type read or modified 
 *
 *********************************************************************
 */

void MethodPpm::add_argument(std::string name) throw()
{
}

/**
 *********************************************************************
 *
 * @return        There is no return value
 *
 * Apply the PPM method
 *
 *********************************************************************
 */

void MethodPpm::apply() throw()
{
}
    
