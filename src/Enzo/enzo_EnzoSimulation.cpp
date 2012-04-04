// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulation user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulation::EnzoSimulation
(
 const char *   parameter_file,
#ifdef CONFIG_USE_CHARM
 int n,
 CProxy_BlockReduce proxy_block_reduce
#else
 GroupProcess * group_process
#endif
 ) throw()
  : Simulation
#ifdef CONFIG_USE_CHARM
    (parameter_file, n, proxy_block_reduce)
#else
    (parameter_file, group_process)
#endif
{
  problem_ = new EnzoProblem;
}

#ifdef CONFIG_USE_CHARM
/// Initialize an empty EnzoSimulation
EnzoSimulation::EnzoSimulation() 
{
};

/// Initialize a migrated EnzoSimulation
EnzoSimulation::EnzoSimulation (CkMigrateMessage *m) 
{
};

//==================================================

#endif

//----------------------------------------------------------------------

EnzoSimulation::~EnzoSimulation() throw()
{
  delete problem_; problem_ = 0;
}

//----------------------------------------------------------------------

void EnzoSimulation::initialize() throw()
{
  // Initialize Cello Simulation base class
  Simulation::initialize();

  // Initialize Enzo-P namespace variables
  enzo::initialize(parameters_,field_descr());

}

//----------------------------------------------------------------------

void EnzoSimulation::finalize() throw()
{
}

//----------------------------------------------------------------------

const Factory * EnzoSimulation::factory() const throw()
{ 
  if (factory_ == NULL) factory_ = new EnzoFactory;
  return factory_;
}


