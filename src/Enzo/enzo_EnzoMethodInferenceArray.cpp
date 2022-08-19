// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodInferenceArray.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-07-12
/// @brief    Implements the EnzoMethodInferenceArray class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodInferenceArray::EnzoMethodInferenceArray
(int level,
 const int root_size[3],
 const int array_dims[3],
 const int array_size[3],
 const int array_ghosts[3],
 std::string field_group)
  : Method(),
    level_(level),
    array_dims_(),
    array_size_(),
    array_ghosts_(),
    field_group_(field_group),
    is_sync_child_(-1),
    is_sync_parent_(-1)
{
  for (int i=0; i<3; i++) {
    array_dims_[i]   = array_dims[i];
    array_size_[i]   = array_size[i];
    array_ghosts_[i] = array_ghosts[i];
  }

  // Compute ratio of cell sizes in level 0 to level "level_"
  int lfactor = 1;
  for (int i=0; i<level_; i++) lfactor*=2;

  // Consistency check: root_size must be array_dims*array_size/(2^level_)

  for (int i=0; i<cello::rank(); i++) {
  ASSERT5("EnzoMethodInferenceArray()",
         "axis-%d parameters mismatch: Root dims [%d] != (array dims [%d]* array size [%d]) / 2^(level [%d])",
          i,root_size[i], array_dims_[i], array_size[i], level_,
          (root_size[i]*lfactor == (array_dims_[i]*array_size[i])));
  }

  cello::define_field ("density");

  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("density");

  // Create Level Arary
  if (CkMyPe() == 0) {

    CProxy_MappingArray array_map = CProxy_MappingArray::ckNew
      (array_dims_[0],array_dims_[1],array_dims_[2]);

    CkArrayOptions opts(array_dims_[0],array_dims_[1],array_dims_[2]);
    opts.setMap(array_map);

    proxy_level_array = CProxy_EnzoLevelArray::ckNew(opts);

    proxy_enzo_simulation.p_set_level_array(proxy_level_array);

  }
  // Create and initialize Block synchronization counters
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  is_sync_child_ = scalar_descr_sync->new_value
    ("method_inference_array_sync_child");
  is_sync_parent_ = scalar_descr_sync->new_value
    ("method_inference_array_sync_parent");
}

//----------------------------------------------------------------------
void EnzoSimulation::p_set_level_array (CProxy_EnzoLevelArray proxy)
{
  proxy_level_array = proxy;
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_level_array_request_data(int ix,int iy,int iz,Index index)
{
  EnzoMethodInferenceArray * method =
    static_cast<EnzoMethodInferenceArray*> (this->method());
  if (is_leaf()) {
    // Pack and send data to inference array (ix,iy,iz)
  } else {
    // Forward request to /overlapping/ child Blocks
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_level_array_send_data(int ix,int iy,int iz)
{
  EnzoMethodInferenceArray * method =
    static_cast<EnzoMethodInferenceArray*> (this->method());

  bool have_all_data = false;
  
  // if leaf block
  if (is_leaf()) {
    // pack data into array
    have_all_data = true;
  } else {
    // accumulate received data
    // sync count /expected/ children: set have_all_data = true if ready
  }

  // If arrays are complete
  if (have_all_data) {
    // if root-level
    if (level() == 0) {
      // send data to inference array element (ix,iy,iz)
    } else {
      // send restricted data to parent Block
    }
  }
}

//----------------------------------------------------------------------

void EnzoBlock:: p_method_level_array_done(int ix,int iy,int iz)
{
  EnzoMethodInferenceArray * method =
    static_cast<EnzoMethodInferenceArray*> (this->method());

  // for each Block overlapping this array
  //    call p_method_level_array_done() on overlapping children
  //    increment Block sync_level_arrays
  //    (initialized with #overlapping inference arrays)
  //    (including non-leaf)
  //    if (sync) compute_done()
}

//----------------------------------------------------------------------

void EnzoMethodInferenceArray::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | level_;
  PUParray(p,array_dims_,3);
  PUParray(p,array_size_,3);
  PUParray(p,array_ghosts_,3);
  p | field_group_;
}

//----------------------------------------------------------------------

void EnzoMethodInferenceArray::compute ( Block * block) throw()
{
  if (block->level() == 0) {
    int n3[3];
    cello::hierarchy()->root_size(n3,n3+1,n3+2);
    int na3[3] = {array_dims_[0],  array_dims_[1],  array_dims_[2]};
    int ga3[3] = {array_ghosts_[0],array_ghosts_[1],array_ghosts_[2]};
    int ib3[3];
    int nb3[3];
    cello::hierarchy()->root_blocks(nb3,nb3+1,nb3+2);
    block->index().array(ib3,ib3+1,ib3+2);

    WARNING("EnzoMethodInferenceArray::compute()",
            "Assuming block ghost depth = 0");

    int ia3_lower[3] = {0,0,0};
    int ia3_upper[3] = {1,1,1};
    intersecting_level_arrays_
      (ia3_lower,ia3_upper, ib3,n3, level_,na3,ga3, nb3);
  }

  // (temporary: compute_done() not called until inference completed)
  block->compute_done();
}

//======================================================================

/* readonly */ int array_size[3];

void EnzoMethodInferenceArray::intersecting_root_blocks_
(
 int ib3_lower[3], // OUTPUT lower-indicies of intersecting root-blocks
 int ib3_upper[3],  // OUTPUT upper-indicies of intersecting root-blocks
 const int ia3[3], // Input: index of intersecting level array
 const int n3[3],  // Input: size of root grid
 const int level,  // Input: level of level array
 const int na3[3], // Input: dimensions of level array (elements)
 const int ga3[3], // Input: ghost depth of inference arrays
 const int nb3[3] // Input: dimensions of root-level mesh
 ) const
{
  CkPrintf ("DEBUG_LEVEL_ARRAY root grid size  %d %d %d\n",n3[0],n3[1],n3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY level array level %d\n",level);
  CkPrintf ("DEBUG_LEVEL_ARRAY array dims na3 %d %d %d\n",na3[0],na3[1],na3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array ghost ga3 %d %d %d\n",ga3[0],ga3[1],ga3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array index ia3 %d %d %d\n",ia3[0],ia3[1],ia3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY block dims nb3 %d %d %d\n",nb3[0],nb3[1],nb3[2]);

  int lfactor = 1; for (int i=0; i<level; i++) lfactor*=2;
  int ma3[3],mb3[3];
  for (int i=0; i<cello::rank(); i++) {
    ma3[i] = n3[i]/na3[i]; // number of root-level cells per inference array
    mb3[i] = n3[i]/nb3[i]; // size of root-level blocks
    //    ia3_lower[i] = std::floor(float((ib3[i]*mb3[i]-ga3[i]-ma3[i])/ma3[i]))+1;
    ib3_lower[i] = std::floor(float((ia3[i]*ma3[i]-ga3[i]))/mb3[i]);
    ib3_upper[i] = std::ceil (float((ia3[i]*ma3[i]+ga3[i]+ma3[i]-1)/mb3[i]))+1;
  }

  CkPrintf ("DEBUG_LEVEL_ARRAY block lower index ib3_lower %d %d %d\n",ib3_lower[0],ib3_lower[1],ib3_lower[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY block upper index ib3_upper %d %d %d\n",ib3_upper[0],ib3_upper[1],ib3_upper[2]);
}

//----------------------------------------------------------------------

void EnzoMethodInferenceArray::intersecting_level_arrays_
(
 int ia3_lower[3], // OUTPUT lower-indicies of intersecting level arrays
 int ia3_upper[3], // OUTPUT upper-indicies of intersecting level arrays
 const int ib3[3], // Input: index of intersecting root-level block
 const int n3[3],  // Input: size of root grid
 const int level,  // Input: level of level array
 const int na3[3], // Input: dimensions of level array (elements)
 const int ga3[3],  // Input: ghost depth of inference arrays
 const int nb3[3] // Input: dimensions of root-level mesh
 ) const
{
  // compute level factor
  int lfactor = 1; for (int i=0; i<level; i++) lfactor*=2;
  int ma3[3],mb3[3];
  for (int i=0; i<cello::rank(); i++) {
    ma3[i] = n3[i]/na3[i]; // number of root-level cells per inference array
    mb3[i] = n3[i]/nb3[i]; // size of root-level blocks
    //    ia3_lower[i] = std::floor(float((ib3[i]*mb3[i]-ga3[i]-ma3[i])/ma3[i]))+1;
    ia3_lower[i] = std::floor(float((ib3[i]*mb3[i]-ga3[i]))/ma3[i]);
    ia3_upper[i] = std::ceil (float((ib3[i]*mb3[i]+ga3[i]+mb3[i]-1)/ma3[i]))+1;
  }
  CkPrintf ("DEBUG_LEVEL_ARRAY root grid size  %d %d %d\n",n3[0],n3[1],n3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY level array level %d\n",level);
  CkPrintf ("DEBUG_LEVEL_ARRAY block size  mb3 %d %d %d\n",mb3[0],mb3[1],mb3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY block dims  nb3 %d %d %d\n",nb3[0],nb3[1],nb3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY block index ib3 %d %d %d\n",ib3[0],ib3[1],ib3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array size  ma3 %d %d %d\n",na3[0],na3[1],na3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array dims  na3 %d %d %d\n",na3[0],na3[1],na3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array ghost ga3 %d %d %d\n",ga3[0],ga3[1],ga3[2]);

  CkPrintf ("DEBUG_LEVEL_ARRAY array lower index ia3_lower %d %d %d\n",ia3_lower[0],ia3_lower[1],ia3_lower[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array upper index ia3_upper %d %d %d\n",ia3_upper[0],ia3_upper[1],ia3_upper[2]);
}

