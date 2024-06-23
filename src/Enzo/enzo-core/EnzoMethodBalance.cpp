// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodBalance.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     <2022-12-09 Fri>
/// @brief    Implements the EnzoMethodBalance class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodBalance::EnzoMethodBalance()
  : Method()
{
  cello::define_field("density");

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("density");
}

//----------------------------------------------------------------------

void EnzoMethodBalance::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
}

//----------------------------------------------------------------------

void EnzoMethodBalance::compute ( Block * block) throw()
{
  // Output that we're calling the load balancer
  if (block->index().is_root())
    cello::monitor()->print("Method EnzoMethodBalance", "entering Cello load-balancer");

  long long index, count;
  block->get_order (&index,&count);

  // Use long long to prevent overflow on large runs
  int ip_next = (long long) CkNumPes()*index/count;

  block->set_ip_next(ip_next);

  int count_local = 0;
  if (ip_next != CkMyPe()) {
    count_local = 1;
  }

  CkCallback callback
    (CkIndex_EnzoSimulation::r_method_balance_count(nullptr), 0,
     proxy_enzo_simulation);

  block->contribute(sizeof(int), &count_local,
                    CkReduction::sum_int, callback);

}

void EnzoSimulation::r_method_balance_count(CkReductionMsg * msg)
{
  // Include self-call of balance check to prevent hanging when no
  // blocks migrate
  int * count_total = (int * )msg->getData();
  sync_method_balance_.set_stop(*count_total + 1);

  if (CkMyPe() == 0) {
    cello::monitor()->print("Method EnzoMethodBalance", "migrating %d blocks",*count_total);
  }

  // Initiate migration
  enzo::block_array().p_method_balance_migrate();
  p_method_balance_check();
}

void EnzoBlock::p_method_balance_migrate()
{
  static_cast<EnzoMethodBalance*> (method())->do_migrate(this);
}

void EnzoMethodBalance::do_migrate(EnzoBlock * enzo_block)
{
  int ip_next = enzo_block->ip_next();
  enzo_block->set_ip_next(-1);
  if (ip_next != CkMyPe()) {
    enzo_block->migrateMe(ip_next);
  }
}

void EnzoSimulation::p_method_balance_check()
{
  if (sync_method_balance_.next()) {
    sync_method_balance_.reset();
    enzo::block_array().doneInserting();
    enzo::block_array().p_method_balance_done();
    cello::monitor()->print("Method EnzoMethodBalance", "done migrating");
  }
}

void EnzoBlock::p_method_balance_done()
{
  static_cast<EnzoMethodBalance*> (method())->done(this);
}

void EnzoMethodBalance::done(EnzoBlock * enzo_block)
{
  enzo_block->compute_done();
}

