// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodInference.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-07-12
/// @brief    Implements the EnzoMethodInference class

#include "cello.hpp"
#include "enzo.hpp"

#define TRACE_INFER

#ifdef TRACE_INFER
#  define TRACE_INFER_BLOCK(B,MSG) \
  CkPrintf ("TRACE_INFER_BLOCK %s %s\n", \
            B->name().c_str(),std::string(MSG).c_str());
#  define TRACE_INFER_SIMULATION(MSG) \
  CkPrintf ("TRACE_INFER_SIMULATION %d %s\n",   \
            CkMyPe(),std::string(MSG).c_str());
#else
#  define TRACE_INFER_BLOCK(B,MSG) /* ... */
#  define TRACE_INFER_SIMULATION(SIM,MSG) /* ... */
#endif


//----------------------------------------------------------------------

EnzoMethodInference::EnzoMethodInference
(int level_base,
 int level_array,
 int level_infer,
 const int array_dims[3],
 const int array_size[3],
 std::string field_group,
 int index_refine,
 int num_refine)
  : Method(),
    level_base_(level_base),
    level_array_(level_array),
    level_infer_(level_infer),
    m3_level_(),
    m3_infer_(),
    field_group_(field_group),
    is_sync_child_(-1),
    is_sync_parent_(-1),
    is_mask_(-1),
    index_refine_(index_refine),
    num_refine_(num_refine)
{
  for (int i=0; i<3; i++) {
    m3_level_[i] = array_dims[i];
    m3_infer_[i] = array_size[i];
  }

  // Verify level_array and array_dims match
  int nd3[3] = {
    cello::config()->mesh_root_blocks[0],
    cello::config()->mesh_root_blocks[1],
    cello::config()->mesh_root_blocks[2] };
  int md3[3] = {
    cello::config()->mesh_root_size[0],
    cello::config()->mesh_root_size[1],
    cello::config()->mesh_root_size[2] };

  
  for (int axis=0; axis<cello::rank(); axis++) {
    ASSERT2("EnzoMethodInference()",
            "level-array level %d mismatch with level array size %d",
            level_array,m3_level_[axis],pow(2,level_array)*nd3[axis]==m3_level_[axis]);
  }

  for (int axis=0; axis<cello::rank(); axis++) {
    int nb = md3[axis]/nd3[axis]; // compute block size along axis
    ASSERT5("EnzoMethodInference()",
            "2**(%d - %d)*(%d/%d) "
            "mismatch with inference array size %d",
            level_infer,level_array,md3[axis],nd3[axis],m3_infer_[axis],
            pow(2,(level_infer-level_array))*nb==m3_infer_[axis]);
  }

  cello::define_field ("density");

  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("density");

  // Create Level Arary
  if (CkMyPe() == 0) {

    //    CProxy_MappingArray array_map = CProxy_MappingArray::ckNew
    //      (m3_level_[0],m3_level_[1],m3_level_[2]);

    //    CkArrayOptions opts(m3_level_[0],m3_level_[1],m3_level_[2]);
    //    opts.setMap(array_map);

    CkPrintf ("Calling EnzoLevelArray::ckNew()\n");
    //    proxy_level_array = CProxy_EnzoLevelArray::ckNew(opts);
    proxy_level_array = CProxy_EnzoLevelArray::ckNew();
    CkPrintf ("Done\n");

    proxy_enzo_simulation.p_set_level_array(proxy_level_array);
    Index3 index0(0,0,0);
    CkPrintf ("Inserting...\n");
    proxy_level_array[index0].insert();
    proxy_level_array.doneInserting();

  }
  // Create and initialize Block synchronization counters
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  is_sync_child_ = scalar_descr_sync->new_value
    ("method_inference_sync_child");
  is_sync_parent_ = scalar_descr_sync->new_value
    ("method_inference_sync_parent");
  is_mask_ = scalar_descr_void->new_value
    ("method_inference_mask");
}

//----------------------------------------------------------------------
void EnzoSimulation::p_set_level_array (CProxy_EnzoLevelArray proxy)
{
  proxy_level_array = proxy;
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_inference_request_data(int ix,int iy,int iz,Index index)
{
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());
  if (is_leaf()) {
    // Pack and send data to inference array (ix,iy,iz)
  } else {
    // Forward request to /overlapping/ child Blocks
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_inference_send_data(int ix,int iy,int iz)
{
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

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

void EnzoBlock::p_method_inference_done(int ix,int iy,int iz)
{
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

  // for each Block overlapping this array
  //    call p_method_inference_done() on overlapping children
  //    increment Block sync_inferences
  //    (initialized with #overlapping inference arrays)
  //    (including non-leaf)
  //    if (sync) compute_done()
}

//----------------------------------------------------------------------

void EnzoMethodInference::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | level_base_;
  PUParray(p,m3_level_,3);
  PUParray(p,m3_infer_,3);
  p | field_group_;
}

//----------------------------------------------------------------------

void EnzoMethodInference::compute ( Block * block) throw()
{
  TRACE_INFER_BLOCK(block,"compute()");
  const int level = block->level();

  // Inialize block scalars
  const int nc = cello::num_children();
  sync_child_ (block).set_stop(nc+1);
  sync_parent_(block).set_stop(1+1);
  *scalar_mask_(block) = nullptr;

  // Negative levels are not involved in EnzoMethodInference
  if (block->level() < 0) {
    TRACE_INFER_BLOCK(block,"compute_done()");
    block->compute_done();
    return;
  }

  if (block->is_leaf()) {

    //    EnzoMsgInferCreate * msg = new EnzoMsgInferCreate;

    // Apply inference array creation criteria
    int n = 0;
    char * mask = nullptr;
    int count = apply_criteria_(block,&n, &mask);

    forward_create_array_(block,count,n,mask);
  }

  // (temporary: compute_done() not called until inference completed)
  TRACE_INFER_BLOCK(block,"compute_done()");
  block->compute_done();
}

//----------------------------------------------------------------------

int EnzoMethodInference::apply_criteria_
(Block * block, int * n, char ** mask)
{
  TRACE_INFER_BLOCK(block,"apply_criteria()");
  // Apply criteria whether to create the inference array under
  // any cells. Return mask of inference arrays to create and
  // the total count
  int count = 0;
  int mask_size = 0;
  if (block->level() >= level_base_) {
    TRACE_INFER_BLOCK(block,"apply_criteria() TRUE");
    // TEMPORARY: create all underlying inference arrays
    const int nc = cello::num_children();
    mask_size = int(pow(nc,(level_array_ - block->level())));
    mask_size = std::max(mask_size,1);
    *n = mask_size;
    *mask = new char[mask_size];
    std::fill_n(*mask,*n,1);
    count += mask_size;
  }
  return count;
}

//----------------------------------------------------------------------

void EnzoMethodInference::forward_create_array_
(Block * block, int count, int n, char * mask)
{
  // Called by (all) leaf blocks
  ASSERT ("EnzoMethodInference::forward_create_array_",
          "Block %s is not a leaf",
          block->is_leaf());

  char buffer[80];
  sprintf (buffer,"forward_create_array_(): leaf %d count %d",
           block->is_leaf()?1:0,count);
  TRACE_INFER_BLOCK(block,buffer);

  int level = block->level();

  // Create level arrays if this is the base level
  if (level == level_base_) {
    create_level_arrays_(block,count,n,mask);
  }

  // Forward count and mask to parent to merge

  if (level == 0) {
    // Root level: synchronize in root simulation object
    proxy_enzo_simulation[0].p_infer_count_arrays(count);
  } else { // level > 0
    int ic3[3] = {0,0,0};
    Index index_parent = block->index().index_parent();
    block->index().child(level,ic3,ic3+1,ic3+2);
    if (level > level_base_) {
      enzo::block_array()[index_parent].
        p_method_infer_merge_masks (count,n,mask,ic3);
      // Count level arrays by beginning reduction to root Simulation
    } else { // 0 < level <= level_bask
      // Forward count to root to count
      enzo::block_array()[index_parent].
        p_method_infer_count_arrays (count);
    }
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_infer_count_arrays(int count)
{
  char buffer[80];
  sprintf(buffer,"p_infer_count_arrays %d",count);
  
  TRACE_INFER_SIMULATION(buffer);
}

//----------------------------------------------------------------------
void EnzoBlock::p_method_infer_merge_masks
(int count, int n, char *mask, int ic3[3])
{
  TRACE_INFER_BLOCK(this,"p_method_infer_merge_masks");

  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());
  CkPrintf ("Method = %p\n",(void *)method);
  method->merge_masks(this,count,n,mask,ic3);
}

void EnzoMethodInference::merge_masks
(Block * block, int count, int n, char *mask, int ic3[3])
{
  char buffer[80];
  sprintf (buffer,"merge_masks() count %d n %d\n",count,n);
  TRACE_INFER_BLOCK (block,buffer);
  if (sync_child_(block).next()) {
    TRACE_INFER_BLOCK (block,"Sync release");
  }
}


//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_count_arrays (int count)
{
  TRACE_INFER_BLOCK(this,"p_method_infer_count_arrays");

  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

  method->count_arrays(this,count);
}

void EnzoMethodInference::count_arrays (Block * block, int count)
{
  char buffer[80];
  sprintf (buffer,"count_arrays() count %d\n",count);
  TRACE_INFER_BLOCK (block,buffer);
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//======================================================================

void EnzoMethodInference::create_level_arrays_
(Block * block, int count, int n, char * mask)
{
  TRACE_INFER_BLOCK(block,"create_level_arrays_()");
}
//----------------------------------------------------------------------
void EnzoMethodInference::intersecting_root_blocks_
(
 int ib3_lower[3], // OUTPUT lower-indicies of intersecting root-blocks
 int ib3_upper[3],  // OUTPUT upper-indicies of intersecting root-blocks
 const int ia3[3], // Input: index of intersecting level array
 const int n3[3],  // Input: size of root grid
 const int level,  // Input: level of level array
 const int na3[3], // Input: dimensions of level array (elements)
 const int nb3[3] // Input: dimensions of root-level mesh
 ) const
{
  CkPrintf ("DEBUG_LEVEL_ARRAY root grid size  %d %d %d\n",n3[0],n3[1],n3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY level array level %d\n",level);
  CkPrintf ("DEBUG_LEVEL_ARRAY array dims na3 %d %d %d\n",na3[0],na3[1],na3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array index ia3 %d %d %d\n",ia3[0],ia3[1],ia3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY block dims nb3 %d %d %d\n",nb3[0],nb3[1],nb3[2]);

  int ma3[3],mb3[3];
  for (int i=0; i<cello::rank(); i++) {
    ma3[i] = n3[i]/na3[i]; // number of root-level cells per inference array
    mb3[i] = n3[i]/nb3[i]; // size of root-level blocks
    //    ia3_lower[i] = std::floor(float((ib3[i]*mb3[i]-ma3[i])/ma3[i]))+1;
    ib3_lower[i] = std::floor(float((ia3[i]*ma3[i]))/mb3[i]);
    ib3_upper[i] = std::ceil (float((ia3[i]*ma3[i]+ma3[i]-1)/mb3[i]))+1;
  }

  CkPrintf ("DEBUG_LEVEL_ARRAY block lower index ib3_lower %d %d %d\n",ib3_lower[0],ib3_lower[1],ib3_lower[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY block upper index ib3_upper %d %d %d\n",ib3_upper[0],ib3_upper[1],ib3_upper[2]);
}

//----------------------------------------------------------------------

void EnzoMethodInference::intersecting_level_arrays_
(
 int ia3_lower[3], // OUTPUT lower-indicies of intersecting level arrays
 int ia3_upper[3], // OUTPUT upper-indicies of intersecting level arrays
 const int ib3[3], // Input: index of intersecting root-level block
 const int n3[3],  // Input: size of root grid
 const int level,  // Input: level of level array
 const int na3[3], // Input: dimensions of level array (elements)
 const int nb3[3] // Input: dimensions of root-level mesh
 ) const
{
  // compute level factor
  int ma3[3],mb3[3];
  for (int i=0; i<cello::rank(); i++) {
    ma3[i] = n3[i]/na3[i]; // number of root-level cells per inference array
    mb3[i] = n3[i]/nb3[i]; // size of root-level blocks
    //    ia3_lower[i] = std::floor(float((ib3[i]*mb3[i]-ma3[i])/ma3[i]))+1;
    ia3_lower[i] = std::floor(float((ib3[i]*mb3[i]))/ma3[i]);
    ia3_upper[i] = std::ceil (float((ib3[i]*mb3[i]+mb3[i]-1)/ma3[i]))+1;
  }
  CkPrintf ("DEBUG_LEVEL_ARRAY root grid size  %d %d %d\n",n3[0],n3[1],n3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY level array level %d\n",level);
  CkPrintf ("DEBUG_LEVEL_ARRAY block size  mb3 %d %d %d\n",mb3[0],mb3[1],mb3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY block dims  nb3 %d %d %d\n",nb3[0],nb3[1],nb3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY block index ib3 %d %d %d\n",ib3[0],ib3[1],ib3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array size  ma3 %d %d %d\n",na3[0],na3[1],na3[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array dims  na3 %d %d %d\n",na3[0],na3[1],na3[2]);

  CkPrintf ("DEBUG_LEVEL_ARRAY array lower index ia3_lower %d %d %d\n",ia3_lower[0],ia3_lower[1],ia3_lower[2]);
  CkPrintf ("DEBUG_LEVEL_ARRAY array upper index ia3_upper %d %d %d\n",ia3_upper[0],ia3_upper[1],ia3_upper[2]);
}

