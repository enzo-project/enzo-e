// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBfieldMethod.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Tues April 20 2021
/// @brief    [\ref Enzo] Implementation of the BfieldMethod abstract base
///           class.

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoBfieldMethod::EnzoBfieldMethod(int num_partial_timesteps)
  : num_partial_timesteps_(num_partial_timesteps),
    partial_timestep_index_(-1),
    target_block_(nullptr)
{
  ASSERT("EnzoConstrainedTransport", "num_partial_timesteps must be positive",
         num_partial_timesteps_ > 0);

  if (num_partial_timesteps_ != 2){
    ERROR("EnzoConstrainedTransport",
          "This machinery hasn't been tested for cases when "
          "num_partial_timesteps!=2.");
  }
}

//----------------------------------------------------------------------

void EnzoBfieldMethod::pup (PUP::er &p)
{
  PUP::able::pup(p);

  p|num_partial_timesteps_;
  if (p.isUnpacking()){
    target_block_ = nullptr;
    // set partial_timestep_index_ to -1 to indicate need for reallocating
    // scratch space
    partial_timestep_index_ = -1;
  } else if (target_block_ != nullptr){
    ERROR("EnzoBfieldMethod::pup",
          "target_block_ is expected to be a nullptr while packing up data");
  }
}

//----------------------------------------------------------------------

void EnzoBfieldMethod::register_target_block(Block *block) noexcept
{
  if (target_block_ != nullptr){
    ERROR("EnzoConstrainedTransport::register_target_block",
          "A target block is already specified");
  }
  target_block_ = block;

  ASSERT("EnzoBfieldMethod::register_target_block", "Invalid state.",
         partial_timestep_index_ <= 0);
  bool first_initialization = (partial_timestep_index_ < 0);

  // Pre-load any method specific data from the new target block and initialize
  // scratch space if this is the first initialization
  register_target_block_(block, first_initialization);

  partial_timestep_index_ = 0;
}

//----------------------------------------------------------------------

void EnzoBfieldMethod::increment_partial_timestep() noexcept
{
  require_registered_block_(); // confirm that target_block_ is valid

  if ((partial_timestep_index_ + 1) == num_partial_timesteps_){
    partial_timestep_index_ = 0;
    target_block_ = nullptr;
  } else {
    partial_timestep_index_++;
  }
}
