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

  // Initialize Simulation synchronization counters
  const int num_root_blocks = nd3[0]*nd3[1]*nd3[2];
  enzo::simulation()->set_sync_infer_count(num_root_blocks);
}

//----------------------------------------------------------------------
void EnzoSimulation::p_set_level_array (CProxy_EnzoLevelArray proxy)
{
  proxy_level_array = proxy;
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

  const int num_blocks_base = cello::hierarchy()->num_blocks(level_base_);
  CkPrintf ("DEBUG_INFER %d num_blocks (level_base) %d\n",CkMyPe(),
            cello::simulation()->num_blocks_level(level_base_));

  // Inialize block scalars
  const int nc = cello::num_children();
  sync_child_ (block).set_stop(nc+1);
  sync_parent_(block).set_stop(1+1);
  *scalar_mask_(block) = nullptr;

  int count = 0;
  if (block->is_leaf()) {

    //    EnzoMsgInferCreate * msg = new EnzoMsgInferCreate;

    // Apply inference array creation criteria
    apply_criteria_(block);

    forward_create_array_(block,count);

  } else if (level >= 0) {

    // Self call merge_masks to ensure block has reached this method
    if (0 <= level && level < level_base_) {
      CkPrintf ("TRACE_COUNT_ARRAYS parent %s\n",block->name().c_str());
      count_arrays(block,0);
    } else {
      int ic3[3] = {-1,-1,-1}; // not accessed
      int count,n;
      char * mask = *scalar_mask_(block);
      merge_masks (block, n=0, mask=nullptr, ic3);
    }
  }
}

//----------------------------------------------------------------------

int EnzoMethodInference::apply_criteria_ (Block * block)
{
  TRACE_INFER_BLOCK(block,"apply_criteria()");
  // Apply criteria whether to create the inference array under
  // any cells. Return mask of inference arrays to create and
  // the total count
  int mask_size = 0;
  char ** mask = scalar_mask_(block);
  int count = 0;
  if (block->level() >= level_base_) {
    // allocate mask
    mask_size = mask_allocate_(block);
    // initialize mask
    int value = (block->level() >= level_base_ + 1) ? 1 : 0;
    //int value = 1;
    std::fill_n(*mask,mask_size,value);
    count += value*mask_size;
  }
  CkPrintf ("mask_size %d\n",mask_size);
  CkPrintf ("Block mask %p\n",*scalar_mask_(block));
  CkPrintf ("Block apply count %s %d\n",block->name().c_str(),count);
  //  print_mask_(block);
  return count;
}

//----------------------------------------------------------------------

int EnzoMethodInference::mask_allocate_(Block * block)
{
  // allocate level array mask if not already allocated
  char ** mask = scalar_mask_(block);
  int m = mask_dims_(block->level());
  if (*mask == nullptr) {
    *mask = new char[m];
    std::fill_n(*mask,m,0);
  }
  return m;
}

int EnzoMethodInference::mask_dims_(int level,int *mx, int *my, int *mz) const
{
  const int nc = cello::num_children();
  int mask_size = int(pow(nc,(level_array_ - level)));
  const int rank = cello::rank();
  int len = pow(2,(level_array_ - level));
  if (mx) *mx = (rank >= 1) ? len : 1;
  if (my) *my = (rank >= 2) ? len : 1;
  if (mz) *mz = (rank >= 3) ? len : 1;
  mask_size = std::max(mask_size,1);
  if (mx && my && mz) {
    ASSERT4 ("EnzoMethodInference::mask_dims_",
             "Unexpected mask dims %d %d %d != size %d",
             *mx,*my,*mz,mask_size,
             (*mx)*(*my)*(*mz) == mask_size);
  }
  return mask_size;
  
}

void EnzoMethodInference::forward_create_array_ (Block * block, int count)
{
  // Called by all blocks level >= 0

  char buffer[80];
  sprintf (buffer,"forward_create_array_(): level %d leaf %d count %d",
           block->level(),block->is_leaf()?1:0,count);
  TRACE_INFER_BLOCK(block,buffer);

  int level = block->level();

  // Create level arrays if this is the base level
  if (level == level_base_) {
    create_level_arrays_(block);
  }

  // Forward count and mask to parent to merge

  if (level == 0) {
    // if root level send final counts to root simulation
    proxy_enzo_simulation[0].p_infer_count_arrays(count);
  } else {
    Index index_parent = block->index().index_parent();
    if (level > level_base_) {
      // if level > level_base, accumulate level array masks
      int ic3[3] = {0,0,0};
      block->index().child(level,ic3,ic3+1,ic3+2);
      int n = mask_dims_(block->level());
      char ** mask = scalar_mask_(block);
      // merge masks
      enzo::block_array()[index_parent].
        p_method_infer_merge_masks (n,*mask,ic3);
    } else { // 0 < level <= level_bask
      // if 0 < level <= level_base_, forward count to root
      enzo::block_array()[index_parent].
        p_method_infer_count_arrays (count);
    }
  }
}

//----------------------------------------------------------------------
void EnzoBlock::p_method_infer_merge_masks (int n, char *mask, int ic3[3])
{
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());
  CkPrintf ("Method = %p\n",(void *)method);
  method->merge_masks(this,n,mask,ic3);
}

void EnzoMethodInference::merge_masks
(Block * block, int n_in, char *mask_in, int ic3[3])
{
  char buffer[80];
  sprintf(buffer,"merge_masks() %s",block->name().c_str());
  TRACE_INFER_BLOCK (block,buffer);
  // Merge masks here if child (may be self)
  if (n_in > 0) {
    CkPrintf ("TRACE_INFER merge_masks() %s Merging child %d %d %d\n",
              name().c_str(),ic3[0],ic3[1],ic3[2]);
    int mx,my,mz;
    int n = mask_dims_(block->level(),&mx,&my,&mz);
    // allocate this block's mask if not allocated yet
    char * mask = *scalar_mask_(block);
    if (mask == nullptr) {
      n = mask_allocate_(block);
      mask = *scalar_mask_(block);
    }
    CkPrintf ("TRACE_INFER merge_masks() %s mask_dims %d %d %d n_in %d\n",
              name().c_str(),mx,my,mz,n_in);
    const int mcx = std::max(mx/2,1);
    const int mcy = std::max(my/2,1);
    const int mcz = std::max(mz/2,1);
    for (int icz=0; icz<mcz; icz++) {
      int iz = ic3[2]*mcz+icz;
      for (int icy=0; icy<mcy; icy++) {
        int iy = ic3[1]*mcy+icy;
        for (int icx=0; icx<mcx; icx++) {
          int ix = ic3[0]*mcx+icx;
          int i=ix+mx*(iy+my*iz);
          int ic=icx+mcx*(icy+mcy*icz);
          mask[i] = mask_in[ic]?1:0;
        }
      }
    }
  }

  // Check if done
  if (sync_child_(block).next()) {
    print_mask_(block);
    // Recurse on block forwarding to root
    forward_create_array_ (block,0);
  }
}


//----------------------------------------------------------------------

void EnzoSimulation::p_infer_count_arrays(int count)
{
  char buffer[80];
  sprintf(buffer,"p_infer_count_arrays %d",count);
  TRACE_INFER_SIMULATION(buffer);

  if (sync_infer_count_.next()) {
    TRACE_INFER_SIMULATION("count_arrays DONE");
    enzo::block_array().p_method_infer_create_arrays();
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_count_arrays (int count)
{
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

  method->count_arrays(this,count);
}

void EnzoMethodInference::count_arrays (Block * block, int count)
{
  char buffer[80];
  sprintf (buffer,"count_arrays() count %d\n",count);
  TRACE_INFER_BLOCK (block,buffer);
  if (sync_child_(block).next()) {
    TRACE_INFER_BLOCK (block,"Sync release count_arrays");
    CkPrintf ("Block %s is_leaf %d\n",block->name().c_str(),block->is_leaf());
    // Recurse on block forwarding to root
    forward_create_array_ (block,count);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_create_arrays()
{
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());
  method->create_arrays(this);
}

void EnzoMethodInference::create_arrays (Block * block)
{
  TRACE_INFER_BLOCK(this,"create_arrays");
  if (block->level() == level_base_) {
    int mx,my,mz;
    int n = mask_dims_(block->level(),&mx,&my,&mz);
    // allocate this block's mask if not allocated yet
    char * mask = *scalar_mask_(block);
    int count = 0;
    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
        for (int ix=0; ix<mx; ix++) {
          const int i=ix+mx*(iy+my*iz);
          if (mask[i]) ++count;
        }
      }
    }
    CkPrintf ("TRACE_INFER count = %d\n",count);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_exit()
{
  TRACE_INFER_BLOCK(this,"p_method_infer_exit()");
  compute_done();
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//======================================================================

void EnzoMethodInference::create_level_arrays_ (Block * block)
{
  TRACE_INFER_BLOCK(block,"create_level_arrays_()");
  ASSERT2 ("EnzoMethodInference::create_level_arrays_",
           "Block %s not in expected refinement level %d",
           block->name().c_str(),block->level(),
           block->level() == level_base_);
  int mx,my,mz;
  int m = mask_dims_(block->level(),&mx,&my,&mz);
  char * mask = *scalar_mask_(block);
  if (mask != nullptr) {
    int i3[3];
    block->index().index_level (i3, block->level());
    int nx,ny,nz;
    level_array_dims_(&nx,&ny,&nz);
    for (int kz=0; kz<mz; kz++) {
      for (int ky=0; ky<my; ky++) {
        for (int kx=0; kx<mx; kx++) {
          const int k=kx+mx*(ky + my*kz);
          if (mask[k]) {
            CkPrintf ("TRACE_INFER create_level_array block %d %d %d\n",
                      i3[0]*mx+kx,i3[1]*my+ky,i3[2]*mz+kz);
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodInference::level_array_dims_(int *mx, int *my, int *mz)
{
  // get root blocks in domain
  int bx,by,bz;
  cello::hierarchy()->root_blocks(&bx,&by,&bz);
  // number of base-level blocks per root-block
  const int rank = cello::rank();
  int r = std::pow(2,level_base_);
  // return number of base-level blocks in domain
  (*mx) = (rank >= 1) ? bx*r : 1;
  (*my) = (rank >= 2) ? by*r : 1;
  (*mz) = (rank >= 3) ? bz*r : 1;
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

