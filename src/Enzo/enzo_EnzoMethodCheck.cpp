// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodCheck.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2202-02-12
/// @brief    Implements the EnzoMethodCheck class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm.hpp"
#include "charm_enzo.hpp"
//----------------------------------------------------------------------

EnzoMethodCheck::EnzoMethodCheck ()
  : Method()
{
}

//----------------------------------------------------------------------

void EnzoMethodCheck::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodCheck::compute ( Block * block) throw()
{
  CkPrintf ("TRACE_METHOD_CHECK %s compute()\n",block->name().c_str());

  CkCallback callback(CkIndex_EnzoSimulation::r_method_check_enter(NULL),0,
                      proxy_enzo_simulation);
  int num_files = 10;
  block->contribute(1*sizeof(int), &num_files, CkReduction::max_int, callback);
}

//----------------------------------------------------------------------

void EnzoSimulation::r_method_check_enter(CkReductionMsg *msg)
{
  int num_files = *((int*)msg->getData());
  delete msg;
  CkPrintf ("%d EnzoSimulation::r_method_check_enter() num_files = %d\n",CkMyPe(),num_files);

  // proxy_simulation = proxy_enzo_simulation = CProxy_EnzoSimulation::ckNew
  //   (parameter_file, strlen(parameter_file)+1);

  /// Initialize synchronization counters
  sync_check_writer_created_.set_stop(num_files);
  sync_check_done_.          set_stop(num_files);

  // Create IoEnzoWriter array. Synchronizes by calling
  // EnzoSimulation[0]::p_writer_created() when done
  proxy_io_enzo_writer = CProxy_IoEnzoWriter::ckNew(num_files);
}


//----------------------------------------------------------------------

IoEnzoWriter::IoEnzoWriter() throw ()
  : CBase_IoWriter()
{
  CkPrintf ("%d IoEnzoWriter()\n",CkMyPe());
  std::string ordering = "order_morton";
  proxy_enzo_simulation[0].p_writer_created(10,ordering);
}

//----------------------------------------------------------------------

void EnzoSimulation::p_writer_created(int num_files, std::string ordering)
{
  CkPrintf ("%d EnzoSimulation::p_writer_created()\n",CkMyPe());
  if (sync_check_writer_created_.next()) {
    CkPrintf ("All Created\n");
    enzo::block_array().p_check_write_first(num_files, ordering);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_check_write_first(int num_files, std::string ordering)
{
  ScalarData<int> * scalar_data = data()->scalar_data_int();
  ScalarDescr *     scalar_descr = cello::scalar_descr_int();
  int is_index = scalar_descr->index(ordering+":index");
  int index    = *scalar_data->value(scalar_descr,is_index);
  CkPrintf ("EnzoBlock::p_check_write_first() %s num %d order %s index %d\n",
            name().c_str(),num_files,ordering.c_str(),index);
}
