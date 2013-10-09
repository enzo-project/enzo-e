// See LICENSE_CELLO file for license and copyright information

/// @file     Io_IoLayout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-28
/// @brief    Implementation of the IoLayout class

#include "io.hpp"

//----------------------------------------------------------------------

IoLayout::IoLayout(const Layout * layout) throw ()
  : Io(3,0),
    layout_((Layout*)layout)
{
  TRACE1("layout = %p",layout);
  meta_name_.push_back("layout_process_first");
  meta_name_.push_back("layout_process_count");
  meta_name_.push_back("layout_block_count");

  ASSERT2("IoLayout::IoLayout()",
	 "meta_name.size() [%d] !=  meta_count_ [%d]",
	  meta_name_.size(),meta_count(),
	  meta_name_.size()==meta_count());
}


//----------------------------------------------------------------------

void IoLayout::meta_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{
  Io::meta_value(index,buffer,name,type,nxd,nyd,nzd);

  TRACE1("layout=%p",layout_);
  if (index == 0) {

    TRACE2("process_first_ = %p %d",layout_, layout_->process_first_);
    *buffer = (void *) &layout_->process_first_;
    *type   = scalar_type_int;
    
  } else if (index == 1) {

    TRACE2("process_first_ = %p %d",layout_, layout_->process_count_);
    *buffer = (void *) &layout_->process_count_;
    *type   = scalar_type_int;

  } else if (index == 2) {

    TRACE2("process_first_ = %p %d",layout_,layout_->block_count_[0]);
    *buffer = (void *) layout_->block_count_;
    *type   = scalar_type_int;
    *nxd     = 3;

  }
  
}

//----------------------------------------------------------------------

void IoLayout::data_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
}

//----------------------------------------------------------------------
