// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeCoolingTime.cpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     2019-05-07
/// @brief    Implements the EnzoComputeCoolingTime class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputeCoolingTime::EnzoComputeCoolingTime
()
  : Compute()
{
}

//----------------------------------------------------------------------

void EnzoComputeCoolingTime::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Compute::pup(p);
}

//----------------------------------------------------------------------

void EnzoComputeCoolingTime::compute ( Block * block) throw()
{

  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  enzo_float * ct = field.is_field("cooling_time") ?
                    (enzo_float*) field.values("cooling_time", i_hist_) : NULL;

  if (!ct) {
    ERROR("EnzoComputeCoolingTime::compute()",
          " 'cooling_time' field is not defined as a permanent field");
  }

  compute(block, ct);
}

//---------------------------------------------------------------------

void EnzoComputeCoolingTime::compute ( Block * block, enzo_float * ct) throw()
{

  if (!block->is_leaf()) return;

  compute_(block, ct);
}

//----------------------------------------------------------------------

void EnzoComputeCoolingTime::compute_(Block * block,
                                      enzo_float * ct,
                                      grackle_field_data * grackle_fields /* NULL */
                                    )
{
  const EnzoMethodGrackle* grackle_method = enzo::grackle_method();
  ASSERT("EnzoComputeCoolingTime::compute_()",
         "Grackle must be enabled in order to compute the cooling time",
         grackle_method != nullptr);
  grackle_method->calculate_cooling_time(EnzoFieldAdaptor(block, i_hist_), ct,
                                         0, grackle_fields);
}
