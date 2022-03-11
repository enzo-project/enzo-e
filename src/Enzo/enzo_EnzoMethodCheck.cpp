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

EnzoMethodCheck::EnzoMethodCheck (int num_files, std::string ordering)
  : Method(),
    num_files_(num_files),
    ordering_(ordering)
{
}

//----------------------------------------------------------------------

void EnzoMethodCheck::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | num_files_;
  p | ordering_;

}

//----------------------------------------------------------------------

void EnzoMethodCheck::compute ( Block * block) throw()
{
  CkCallback callback(CkIndex_EnzoSimulation::r_method_check_enter(NULL),0,
                      proxy_enzo_simulation);
  block->contribute(callback);
}

//----------------------------------------------------------------------

void EnzoSimulation::r_method_check_enter(CkReductionMsg *msg)
{
  // [ Called on ip=0 only ]
  
  delete msg;
  // proxy_simulation = proxy_enzo_simulation = CProxy_EnzoSimulation::ckNew
  //   (parameter_file, strlen(parameter_file)+1);

  check_num_files_ =        enzo::config()->method_check_num_files;
  check_ordering_  = enzo::config()->method_check_ordering;

  /// Initialize synchronization counters
  sync_check_writer_created_.set_stop(check_num_files_);
  sync_check_done_.          set_stop(check_num_files_);

  // Create IoEnzoWriter array. Synchronizes by calling
  // EnzoSimulation[0]::p_writer_created() when done

  proxy_io_enzo_writer = CProxy_IoEnzoWriter::ckNew
    (check_num_files_, check_ordering_, check_num_files_);
  proxy_io_enzo_writer.doneInserting();
}

//----------------------------------------------------------------------

IoEnzoWriter::IoEnzoWriter(int num_files, std::string ordering) throw ()
  : CBase_IoEnzoWriter(),
    num_files_(num_files),
    ordering_(ordering)
{
  proxy_enzo_simulation[0].p_writer_created(thisArrayID);
}

//----------------------------------------------------------------------

void EnzoSimulation::p_writer_created (CProxy_IoEnzoWriter proxy)
{
  if (sync_check_writer_created_.next()) {
    proxy_enzo_simulation.p_check_continue (proxy);
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_check_continue (CProxy_IoEnzoWriter proxy)
{
  // Save proxy created on pe 0 on all processes
  proxy_io_enzo_writer = proxy;

  CkCallback callback(CkIndex_EnzoSimulation::r_check_write(NULL),0,
                      proxy_enzo_simulation);
  contribute(callback);
}

//----------------------------------------------------------------------

void EnzoSimulation::r_check_write(CkReductionMsg *msg)
{
  delete msg;
  enzo::block_array().p_check_write_first(check_num_files_, check_ordering_);
}
//----------------------------------------------------------------------

void EnzoBlock::p_check_write_first(int num_files, std::string ordering)
{
  Index index_this, index_next;
  std::string name_this, name_next;
  int index_block, index_file;
  bool is_first, is_last;

  check_get_parameters_
    (num_files,ordering,
     index_this,index_next,
     name_this, name_next,
     index_block,index_file,
     is_first, is_last);

  if (is_first) {
    proxy_io_enzo_writer[index_file].p_write
      (index_file, name_this,name_next, index_this,index_next,
       index_block,is_last);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_check_write_next(int num_files, std::string ordering)
{
  Index index_this, index_next;
  std::string name_this, name_next;
  int index_block, index_file;
  bool is_first, is_last;
  
  check_get_parameters_
    (num_files,ordering,
     index_this,index_next,
     name_this, name_next,
     index_block,index_file,
     is_first, is_last);

  proxy_io_enzo_writer[index_file].p_write
    (index_file,name_this, name_next, index_this,index_next,index_block,is_last);
}

//----------------------------------------------------------------------

void IoEnzoWriter::p_write
(int index_file, std::string name_this, std::string name_next,
 Index index_this, Index index_next,
 int index_block, bool is_last)
{
  if (!is_last) {
    enzo::block_array()[index_next].p_check_write_next(num_files_, ordering_);
  } else {
    proxy_enzo_simulation[0].p_check_done();
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_check_done()
{
  if (sync_check_done_.next()) {
    enzo::block_array().p_check_done();
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_check_done()
{  compute_done(); }

//======================================================================

void EnzoBlock::check_get_parameters_
( int           num_files,
  std::string   ordering,
  Index &       index_this,
  Index &       index_next,
  std::string & name_this,
  std::string & name_next,
  int &         index_block,
  int &         index_file,
  bool &        is_first,
  bool &        is_last
  )
{
  ScalarData<int> *   scalar_data_int    = data()->scalar_data_int();
  ScalarDescr *       scalar_descr_int   = cello::scalar_descr_int();
  ScalarData<Index> * scalar_data_index  = data()->scalar_data_index();
  ScalarDescr *       scalar_descr_index = cello::scalar_descr_index();
  const int is_index = scalar_descr_int->index(ordering+":index");
  const int is_count = scalar_descr_int->index(ordering+":count");
  const int is_next  = scalar_descr_index->index(ordering+":next");
  const int count    = *scalar_data_int->value(scalar_descr_int,is_count);
  const Index next   = *scalar_data_index->value(scalar_descr_index,is_next);

  const int rank = cello::rank();
  Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  int na3[3];
  hierarchy->root_blocks(na3,na3+1,na3+2);

  const int min_level = hierarchy->min_level();

  index_this = this->index();
  index_next = next;

  name_this = this->name();
  name_next = this->name(index_next);

  index_block    = *scalar_data_int->value(scalar_descr_int,is_index);
  index_file = index_block*num_files/count;

  is_first = (index_block == 0) || (index_block*num_files/count != (index_block-1)*num_files/count);
  is_last  = (index_block*num_files/count != (index_block+1)*num_files/count);

}

