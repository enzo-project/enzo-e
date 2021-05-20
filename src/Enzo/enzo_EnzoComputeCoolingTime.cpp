// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeCoolingTime.cpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     2019-05-07
/// @brief    Implements the EnzoComputeCoolingTime class

#ifdef CONFIG_USE_GRACKLE

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
                                      code_units * grackle_units /* NULL */ ,
                                      grackle_field_data * grackle_fields /* NULL */
                                    )
{
  ASSERT("EnzoComputeCoolingTime::compute_()",
         "Grackle must be enabled in order to compute the cooling time",
         enzo::config()->method_grackle_use_grackle );
  const EnzoMethodGrackle* grackle_method = enzo::grackle_method();
  grackle_method->calculate_cooling_time(block, ct, grackle_units,
					 grackle_fields, i_hist_);
}

#endif
