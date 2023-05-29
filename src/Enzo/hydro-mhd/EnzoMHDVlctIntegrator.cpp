// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMHDVlctIntegrator.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sun May 28 2023
/// @brief    [\ref Enzo] Implementation of EnzoMHDVlctIntegrator

#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"

// place this after #include "enzo.hpp"
#include "EnzoMHDVlctIntegrator.hpp"

//----------------------------------------------------------------------

EnzoMHDVlctIntegrator::~EnzoMHDVlctIntegrator()
{
  delete half_dt_recon_;
  delete full_dt_recon_;
  delete riemann_solver_;
  delete integration_quan_updater_;
}

//----------------------------------------------------------------------

