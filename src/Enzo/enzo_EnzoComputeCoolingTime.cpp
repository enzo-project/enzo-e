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
  EnzoBlock * enzo_block = enzo::block(block);

  Field field = enzo_block->data()->field();

  const EnzoConfig * enzo_config = enzo::config();

  ASSERT("EnzoComputeCoolingTime::compute_()",
         "Grackle must be enabled in order to compute the cooling time",
         enzo_config->method_grackle_use_grackle );

  code_units grackle_units_;
  grackle_field_data grackle_fields_;

  // setup grackle units if they are not already provided
  if (!grackle_units){
    grackle_units = &grackle_units_;
    EnzoMethodGrackle::setup_grackle_units(enzo_block, grackle_units, i_hist_);
  }

  // if grackle fields are not provided, define them
  bool delete_grackle_fields = false;
  if (!grackle_fields){
    grackle_fields  = &grackle_fields_;
    EnzoMethodGrackle::setup_grackle_fields(enzo_block, grackle_fields, i_hist_);
    delete_grackle_fields = true;
  }

  // temperature is returned in units of K
  if (calculate_cooling_time(&grackle_units_, &grackle_fields_, ct) == ENZO_FAIL){
    ERROR("EnzoComputeCoolingTime::compute_()",
          "Error in call to Grackle's compute_temperature routine.\n");
  }

  if (delete_grackle_fields){
    EnzoMethodGrackle::delete_grackle_fields(grackle_fields);
  }
  return;
}

#endif
