/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretionComputeBondiHoyle.cpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @author     John Regan (john.regan@mu.ie)
/// @date
/// @brief      Computes accretion rates according to Bondi-Hoyle model.
///             See Krumholz+ 2004, ApJ, 611, 399 for details.
///

#include "cello.hpp"
#include "enzo.hpp"

//------------------------------------------------------------------

EnzoMethodAccretionComputeBondiHoyle::EnzoMethodAccretionComputeBondiHoyle
(double accretion_radius_cells)
  : EnzoMethodAccretionCompute(accretion_radius_cells)
{

}

//-------------------------------------------------------------------

void EnzoMethodAccretionComputeBondiHoyle::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodAccretionCompute::pup(p); // call parent class pup

  return;
}

//--------------------------------------------------------------------

void EnzoMethodAccretionComputeBondiHoyle::compute (Block * block) throw()
{
 
  if (block->is_leaf()){
    this->compute_(block);
  }
  block->compute_done();

  return;
}

//-------------------------------------------------------------------------

void EnzoMethodAccretionComputeBondiHoyle::compute_(Block * block)

{
  return;
}
