/// See LICENSE_CELLO file for license and copyright information

///
///
///
///

#include "cello.hpp"
#include "enzo.hpp"

EnzoMethodFeedback::EnzoMethodFeedback
()
  : Method()
{
  FieldDescr * field_descr = cello::field_descr();
  const EnzoConfig * enzo_config = enzo::config();

  // Initialize default refresh object
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
                             enzo_sync_id_method_feedback);
  refresh(ir)->add_all_fields();

  return;
}

void EnzoMethodFeedback::pup (PUP::er &p)
{
  /// NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  return;
}

void EnzoMethodFeedback::compute (Block * block) throw()
{
  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

//  if (block->is_leaf()){
    // do something here
//  }

  enzo_block->compute_done();

  return;
}

double EnzoMethodFeedback::timestep (Block * block) const throw()
{
  // In general this is not needed, but could imagine putting timestep
  // limiters in situations where, for example, one would want
  // dt < star_lifetime (or something like that), especially if
  // important things happen throughout the star's lifetime.

  return std::numeric_limits<double>::max();
}
