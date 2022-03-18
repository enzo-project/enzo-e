// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodCheck.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2202-02-12
/// @brief    Implements the EnzoMethodCheck class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm.hpp"
#include "charm_enzo.hpp"
#include <iostream>
#include <iomanip>

// #define TRACE

#ifdef TRACE
#   undef TRACE
#   define TRACE(MSG)                                           \
  CkPrintf ("%d TRACE %s\n",CkMyPe(),std::string(MSG).c_str()); \
  fflush(stdout);
#   define TRACE_BLOCK(MSG,BLOCK)                              \
  CkPrintf ("%d TRACE %s %s\n",CkMyPe(),BLOCK->name().c_str(), \
            std::string(MSG).c_str());                         \
  fflush(stdout);
#else
#   define TRACE(MSG)  /* ... */
#   define TRACE_BLOCK(MSG,BLOCK)  /* ... */
#endif

int Simulation::file_counter_ = 0;

//----------------------------------------------------------------------

EnzoMethodCheck::EnzoMethodCheck
(int num_files, std::string ordering, std::vector<std::string> directory)
  : Method(),
    num_files_(num_files),
    ordering_(ordering),
    directory_(directory)
{
  TRACE("[1] EnzoMethodCheck::EnzoMethodCheck()");
  Refresh * refresh = cello::refresh(ir_post_);
  cello::simulation()->refresh_set_name(ir_post_,name());
  refresh->add_field("density");
  // Create IO writer
  if (CkMyPe() == 0) {

    enzo::simulation()->set_sync_check_writer(num_files_);

    proxy_io_enzo_writer = CProxy_IoEnzoWriter::ckNew
      (num_files, ordering, num_files);

    proxy_io_enzo_writer.doneInserting();

    proxy_enzo_simulation.p_set_io_writer(proxy_io_enzo_writer);

  }

}

//----------------------------------------------------------------------

void EnzoSimulation::p_set_io_writer(CProxy_IoEnzoWriter io_enzo_writer)
{ proxy_io_enzo_writer = io_enzo_writer; }

//----------------------------------------------------------------------

void EnzoMethodCheck::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | num_files_;
  p | ordering_;
  p | directory_;

}

//----------------------------------------------------------------------

void EnzoMethodCheck::compute ( Block * block) throw()
{
  TRACE_BLOCK("[2] EnzoMethodCheck::compute()",block);
  CkCallback callback(CkIndex_EnzoSimulation::r_method_check_enter(NULL),0,
                      proxy_enzo_simulation);
  block->contribute(callback);
}

//----------------------------------------------------------------------

void EnzoSimulation::r_method_check_enter(CkReductionMsg *msg)
  // [ Called on ip=0 only ]
{
  TRACE("[3] EnzoSimulation::r_method_check_enter()");
  
  delete msg;

  check_num_files_  = enzo::config()->method_check_num_files;
  check_ordering_   = enzo::config()->method_check_ordering;
  check_directory_  = enzo::config()->method_check_dir;

  /// Initialize synchronization counters
  sync_check_done_.          set_stop(check_num_files_);

  /// Create the directory

  bool already_exists = false;
  
  std::string name_dir = file_create_dir_(check_directory_,already_exists);

  if (already_exists) {
    // Exit checkpoint if directory already exists
    enzo::block_array().p_check_done();
    
  } else {
    // Else start checkpoint

    // Create hierarchy file if root writer

    std::string name_file = name_dir + "/check.file_list";
    std::ofstream stream_file_list (name_file);

    ASSERT1("IoEnzoWriter",
            "Cannot open hierarchy file %s for writing",
            name_file,stream_file_list);

    stream_file_list << std::setfill('0');
    int max_digits = log(check_num_files_-1)/log(10) + 1;
    stream_file_list << check_num_files_ << std::endl;
    for (int i=0; i<check_num_files_; i++) {
      stream_file_list << "block_data-" << std::setw(max_digits) << i << std::endl;
    }

    enzo::block_array().p_check_write_first
      (check_num_files_, check_ordering_, name_dir);
  }
  // Create IoEnzoWriter array. Synchronizes by calling
  // EnzoSimulation[0]::p_writer_created() when done

}

//----------------------------------------------------------------------

IoEnzoWriter::IoEnzoWriter
(int num_files, std::string ordering) throw ()
  : CBase_IoEnzoWriter(),
    num_files_(num_files),
    ordering_(ordering)
{
  TRACE("[4] IoEnzoWriter::IoEnzoWriter()");
}

//----------------------------------------------------------------------

void EnzoBlock::p_check_write_first(int num_files, std::string ordering, std::string name_dir)
{
  TRACE_BLOCK("[8] EnzoBlock::p_check_write_first",this);
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
       index_block,is_first,is_last,name_dir);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_check_write_next(int num_files, std::string ordering)
{
  TRACE_BLOCK("[9] EnzoBlock::p_check_write_next",this);
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
    (index_file,name_this, name_next, index_this,index_next,index_block,
     is_first,is_last,"");
}

//----------------------------------------------------------------------

void IoEnzoWriter::p_write
(int index_file, std::string name_this, std::string name_next,
 Index index_this, Index index_next, int index_block,
 bool is_first, bool is_last, std::string name_dir)
{
  // Write to block list file, opening or closing file as needed

  if (is_first) stream_block_list_ = create_block_list_(name_dir);

  write_block_list_(name_this);
  
  if (is_last) close_block_list_();

  TRACE("[A] IoEnzoWriter::p_write_first");
  if (!is_last) {
    enzo::block_array()[index_next].p_check_write_next(num_files_, ordering_);
  } else {
    proxy_enzo_simulation[0].p_check_done();
  }
}

//----------------------------------------------------------------------

std::ofstream IoEnzoWriter::create_block_list_(std::string name_dir)
{
  // Create hierarchy file if root writer

  char file_name[80];
  sprintf (file_name,"%s/block_data-%d.block_list",
           name_dir.c_str(),thisIndex);

  std::ofstream stream_block_list (file_name);

  ASSERT1("Simulation::create_block_list_",
          "Cannot open block_list file %s for writing",
          file_name,stream_block_list);

  return stream_block_list;
}

//----------------------------------------------------------------------

void IoEnzoWriter::write_block_list_(std::string block_name)
{
  stream_block_list_ << block_name << std::endl;
}

//----------------------------------------------------------------------

void IoEnzoWriter::close_block_list_()
{
}

//----------------------------------------------------------------------

void EnzoSimulation::p_check_done()
{
  TRACE("[B] EnzoSimulation::p_check_done");
  if (sync_check_done_.next()) {
    enzo::block_array().p_check_done();
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_check_done()
{
  TRACE_BLOCK("[C] EnzoBlock::p_check_done",this);
  compute_done();
}

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
  const Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  int na3[3];
  hierarchy->root_blocks(na3,na3+1,na3+2);

  const int min_level = hierarchy->min_level();

  index_this = this->index();
  index_next = next;

  name_this = this->name();
  name_next = this->name(index_next);

  index_block    = *scalar_data_int->value(scalar_descr_int,is_index);
  index_file = index_block*num_files/count;

  const int ib  = index_block;
  const int ibm = index_block - 1;
  const int ibp = index_block + 1;
  const int nb = count;
  const int nf = num_files;
  is_first = (ib  == 0)  || (ib*nf/nb != (ibm)*nf/nb);
  is_last  = (ibp == nb) || (ib*nf/nb != (ibp)*nf/nb);

}

//----------------------------------------------------------------------

std::string Simulation::file_create_dir_
(std::vector<std::string> directory_format, bool & already_exists)
{
  const int counter = Simulation::file_counter_++;
  const int cycle = cello::simulation()->cycle();
  const double time = cello::simulation()->time();

  cello::create_directory
    (&directory_format, counter,cycle,time,already_exists);

  ASSERT("Simulation::file_create_dir_",
         "Directory name must be non-empty",
         directory_format[0] != "");

  std::string name_dir = cello::expand_name
    (&directory_format,counter, cycle, time);

  return name_dir;
}

