// See LICENSE_CELLO file for license and copyright information

/// @file     Io_IoPatch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoPatch class

#include "io.hpp"

//----------------------------------------------------------------------

IoPatch::IoPatch(const Patch * patch) throw ()
  : Io(6,0),
    patch_(patch)
{
  meta_name_.push_back("id");
  meta_name_.push_back("size");
  meta_name_.push_back("offset");
  meta_name_.push_back("blocking");
  meta_name_.push_back("lower");
  meta_name_.push_back("upper");

  ASSERT2("IoPatch::IoPatch()",
	 "meta_name.size() [%d] !=  meta_count_ [%d]",
	  meta_name_.size(),meta_count(),
	  meta_name_.size()==meta_count());
}


//----------------------------------------------------------------------

void IoPatch::meta_value
(int index,
 void ** buffer, std::string * name, enum scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{
  DEBUG0;
  Io::meta_value(index,buffer,name,type,nxd,nyd,nzd);
  DEBUG0;


#ifdef CONFIG_USE_CHARM
  DEBUG0;
  const Patch * patch = ((CProxy_Patch *)patch_)->ckLocal();
  DEBUG0;
#else
  const Patch * patch = patch_;
#endif

  int count = 0;

  DEBUG0;
  DEBUG1 ("patch = %p",patch);
  DEBUG1 ("patch id = %d",patch->id_);
  if (index == count++) {
    *buffer = (void *) &patch->id_;
    *type   = scalar_type_int;
  } else  if (index == count++) {
    *buffer = (void *) patch->size_;
    *type   = scalar_type_int;
    *nxd     = 3;
  } else if (index == count++) {
    *buffer = (void *) patch->offset_;
    *type   = scalar_type_int;
    *nxd     = 3;
  } else if (index == count++) {
    *buffer = (void *) patch->blocking_;
    *type   = scalar_type_int;
    *nxd     = 3;
  } else if (index == count++) {
    *buffer = (void *) patch->lower_;
    *type   = scalar_type_double;
    *nxd     = 3;
  } else if (index == count++) {
    *buffer = (void *) patch->upper_;
    *type   = scalar_type_double;
    *nxd     = 3;
  }
  
}

//----------------------------------------------------------------------

void IoPatch::data_value
(int index,
 void ** buffer, std::string * name, enum scalar_type * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
}

//----------------------------------------------------------------------
