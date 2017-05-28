// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Initial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Problem class

#include "cello.hpp"
#include "main.hpp"
#include "simulation.hpp"
#include "data.hpp"
#include "problem.hpp"

//----------------------------------------------------------------------

extern CProxy_Simulation  proxy_simulation;

//----------------------------------------------------------------------

void Initial::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  PUP::able::pup(p);

  p | cycle_;
  p | time_;

}

//----------------------------------------------------------------------

void Initial::enforce_block
(
 Block            * block,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 const Hierarchy  * hierarchy
 ) throw()
{
  CkPrintf ("Initial::enforce_block %d %g\n",cycle_,time_);
    block->set_cycle(cycle_);
    block->set_time(time_);
}
