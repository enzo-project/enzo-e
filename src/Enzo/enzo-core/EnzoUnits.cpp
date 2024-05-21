// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoUnits.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-07-13
/// @brief    [\ref Enzo] Implementation of the EnzoUnits class

#include "enzo.cpp"

//----------------------------------------------------------------------

void EnzoUnits::pup (PUP::er &p) {

  TRACEPUP;

  Units::pup(p);

  p | cosmology_; // PUP::able
}

//----------------------------------------------------------------------

void EnzoUnits::set_current_time (double time)
{
  if (cosmology_ != nullptr) { cosmology_->set_current_time(time); }
}

//----------------------------------------------------------------------

double EnzoUnits::time() const
{
  return (cosmology_ == nullptr) ?
    Units::time() : cosmology_->time_units();
}

//----------------------------------------------------------------------

double EnzoUnits::mass() const
{
  return (cosmology_ == nullptr) ?
    Units::mass() : cosmology_->mass_units();
}

//----------------------------------------------------------------------

double EnzoUnits::length() const
{
  return (cosmology_ == nullptr) ?
    Units::length() : cosmology_->length_units();
}

//----------------------------------------------------------------------

double EnzoUnits::density() const
{
  return (cosmology_ == nullptr) ?
    Units::density() : cosmology_->density_units();
}

//----------------------------------------------------------------------

double EnzoUnits::velocity() const
{
  return (cosmology_ == nullptr) ?
    Units::velocity() : cosmology_->velocity_units();
}
