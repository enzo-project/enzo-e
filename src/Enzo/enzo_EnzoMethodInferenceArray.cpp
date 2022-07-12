// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodInferenceArray.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-07-12
/// @brief    Implements the EnzoMethodInferenceArray class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodInferenceArray::EnzoMethodInferenceArray
(int level, int array_size, int array_ghosts, std::string field_group)
  : Method(),
    level_(level),
    array_size_(array_size),
    array_ghosts_(array_ghosts),
    field_group_(field_group)
{

  cello::define_field ("density");

  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("density");

}

//----------------------------------------------------------------------

void EnzoMethodInferenceArray::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | level_;
  p | array_size_;
  p | array_ghosts_;
  p | field_group_;
}

//----------------------------------------------------------------------

void EnzoMethodInferenceArray::compute ( Block * block) throw()
{

  if (block->is_leaf()) {
    CkPrintf ("TRACE_INFERENCE_ARRAY %s compute() leaf()\n",
              block->name().c_str());
  }

  CkPrintf ("TRACE_INFERENCE_ARRAY %s compute_done()\n",
            block->name().c_str());

  
  block->compute_done();
}

