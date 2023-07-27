// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodBalance.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     <2022-12-09 Fri>
/// @brief    Implements the EnzoMethodBalance class

#include "cello.hpp"
#include "enzo.hpp"

// #define TRACE_BALANCE

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
#ifdef TRACE_BALANCE
  CkPrintf ("TRACE_SELF_BALANCE %s\n",block->name().c_str());
  CkPrintf ("TRACE_BALANCE 1 compute() %s process %d counter %lld\n",
            block->name().c_str(),CkMyPe(),MsgRefresh::counter[CkMyPe()]);
#endif
  // Output that we're calling the load balancer
  Monitor * monitor = cello::monitor();
  if (block->index().is_root())
    monitor->print("Method", "Calling Cello load-balancer");

  long long index, count;

  block->get_order (&index,&count);
#ifdef TRACE_BALANCE
  CkPrintf ("TRACE_SELF_BALANCE pe %d %s index %d count %d\n",CkMyPe(),block->name().c_str(),
            index,count);
#endif

  // Use long long to prevent overflow on large runs
  int ip_next = (long long) CkNumPes()*index/count;

  block->set_ip_next(ip_next);
#ifdef TRACE_BALANCE
  CkPrintf ("self_balance %lld %lld %d %d\n", count, index,ip_next,CkMyPe());
#endif

  int count_local = 0;
  if (ip_next != CkMyPe()) {
#ifdef TRACE_BALANCE
    CkPrintf ("TRACE_MIGRATE Method Counting %s from %d to %d\n",block->name().c_str(),CkMyPe(),ip_next);
#endif
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
  int * count_total = (int * )msg->getData();
#ifdef TRACE_BALANCE
  CkPrintf ("DEBUG_BALANCE block_count = %d\n",*count_total);
  fflush(stdout);
#endif
  sync_method_balance_.set_stop(*count_total + 1);
  // Initiate migration
  enzo::block_array().p_method_balance_migrate();
  // Include self-call of balance check to prevent hanging of
  // no blocks migrate
  p_method_balance_check();
}

void EnzoBlock::p_method_balance_migrate()
{
  static_cast<EnzoMethodBalance*> (method())->do_migrate(this);
}

void EnzoMethodBalance::do_migrate(EnzoBlock * enzo_block)
{
#ifdef TRACE_BALANCE
  CkPrintf ("TRACE_BALANCE 2 migrate() %s process %d counter %lld\n",
            enzo_block->name().c_str(),CkMyPe(),MsgRefresh::counter[CkMyPe()]);
#endif
  int ip_next = enzo_block->ip_next();
  enzo_block->set_ip_next(-1);
  if (ip_next != CkMyPe()) {
#ifdef TRACE_BALANCE
    CkPrintf ("TRACE_MIGRATE Method Migrating %s from %d to %d\n",
              enzo_block->name().c_str(),CkMyPe(),ip_next);
#endif
    enzo_block->migrateMe(ip_next);
  }
}

void EnzoSimulation::p_method_balance_check()
{
#ifdef TRACE_BALANCE
  CkPrintf ("TRACE_MIGRATE p_method_balance_check()\n");
#endif
#ifdef TRACE_BALANCE
  sync_method_balance_.print("p_method_balance_check()");
  CkPrintf ("TRACE_BALANCE 3 check() process %d counter %lld\n",
            CkMyPe(),MsgRefresh::counter[CkMyPe()]);
  fflush(stdout);
#endif
  if (sync_method_balance_.next()) {
    sync_method_balance_.reset();
    enzo::block_array().doneInserting();
    enzo::block_array().p_method_balance_done();
  }
}

void EnzoBlock::p_method_balance_done()
{
  static_cast<EnzoMethodBalance*> (method())->done(this);
}

void EnzoMethodBalance::done(EnzoBlock * enzo_block)
{
#ifdef TRACE_BALANCE
  CkPrintf ("TRACE_BALANCE 4 done() %s process %d counter %lld\n",
            enzo_block->name().c_str(),CkMyPe(),MsgRefresh::counter[CkMyPe()]);
#endif
  enzo_block->compute_done();
}

