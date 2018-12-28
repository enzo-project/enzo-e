#include "enzo.hpp"
#include "charm_enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialLinearWave::EnzoInitialLinearWave(int cycle, double time,
					     double alpha, double beta,
					     double gamma,
					     int wave_type) throw()
  : Initial (cycle,time),
    alpha_(alpha),
    beta_(beta),
    gamma_(gamma),
    wave_type_(wave_type)
{
}

//----------------------------------------------------------------------

void EnzoInitialLinearWave::enforce_block(Block * block,
					  const Hierarchy * hierarchy) throw()
{
  // Set up the test problem
  // This will only work for the VLCT method
  // (PPML initial conditions are much more complicated)
  // Maybe don't pre-calculate alpha and beta



}
