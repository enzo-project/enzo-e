// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodInference.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-07-12
/// @brief    Implements the EnzoMethodInference class

#include "cello.hpp"
#include "enzo.hpp"

// #define TRACE_SYNC
// #define TRACE_INFER
// #define TRACE_ENTRY
// #define TRACE_ARRAY

#ifdef TRACE_INFER
#undef TRACE_INFER
#  define TRACE_INFER(MSG) \
  CkPrintf ("TRACE_INFER %s\n", \
            std::string(MSG).c_str()); \
  fflush(stdout);
#  define TRACE_INFER_BLOCK(B,MSG) \
  CkPrintf ("TRACE_INFER_BLOCK %s %s\n", \
            B->name().c_str(),std::string(MSG).c_str()); \
  fflush(stdout);
#  define TRACE_INFER_SIMULATION(MSG) \
  CkPrintf ("TRACE_INFER_SIMULATION %d %s\n",   \
            CkMyPe(),std::string(MSG).c_str()); \
  fflush(stdout);
#  define WRITE_ARRAY(A,MX,MY,MZ,NX,NY,NZ,EX,EY,EZ)               \
        {                                                       \
          for (int iz=0; iz<NZ; iz++) {                         \
            for (int iy=0; iy<NY; iy++) {                       \
              for (int ix=0; ix<NX; ix++) {                     \
                auto i = ix + MX*(iy + MY*iz);                  \
                CkPrintf ("%d %d %d %g\n",ix,iy,iz,A[i]);       \
              }                                                 \
            }                                                   \
          }                                                     \
      }
#else
#  define TRACE_INFER(MSG) /* ... */
#  define TRACE_INFER_BLOCK(B,MSG) /* ... */
#  define TRACE_INFER_SIMULATION(MSG) /* ... */
#  define WRITE_ARRAY(A,MX,MY,MZ,NX,NY,NZ,EX,EY,EZ) /* ... */
#endif

#ifdef TRACE_ENTRY
#undef TRACE_ENTRY
#   define TRACE_ENTRY_CALL(ENTRY) CkPrintf ("TRACE_ENTRY Calling %s\n", \
                                             std::string(ENTRY).c_str()); \
  fflush(stdout);
#   define TRACE_ENTRY_START(ENTRY) CkPrintf ("TRACE_ENTRY Starting %s\n", \
                                             std::string(ENTRY).c_str()); \
  fflush(stdout);
#else
#   define TRACE_ENTRY_CALL(ENTRY) /* ... */
#   define TRACE_ENTRY_START(ENTRY) /* ... */
#endif

#ifdef TRACE_ARRAY
#undef TRACE_ARRAY
#   define TRACE_ARRAY(ARRAY) CkPrintf ("TRACE_ARRAY %s\n", \
                                             std::string(ARRAY).c_str()); \
  fflush(stdout);
#   define TRACE_ARRAY_CALL(ARRAY) CkPrintf ("TRACE_ARRAY Calling %s\n", \
                                             std::string(ARRAY).c_str()); \
  fflush(stdout);
#   define TRACE_ARRAY_START(ARRAY) CkPrintf ("TRACE_ARRAY Starting %s\n", \
                                             std::string(ARRAY).c_str()); \
  fflush(stdout);
#else
#   define TRACE_ARRAY(ARRAY) /* ... */
#   define TRACE_ARRAY_CALL(ARRAY) /* ... */
#   define TRACE_ARRAY_START(ARRAY) /* ... */
#endif

#ifdef TRACE_SYNC
#undef TRACE_SYNC
#   define TRACE_SYNC(SYNC,MSG)  SYNC.print(MSG);
#else
#   define TRACE_SYNC(SYNC,MSG) /* ... */
#endif

//----------------------------------------------------------------------

EnzoMethodInference::EnzoMethodInference
(int level_base,
 int level_array,
 int level_infer,
 const int m3_level[3],
 const int m3_infer[3],
 std::string field_group,
 int index_criteria,
 int num_criteria)
  : Method(),
    level_base_(level_base),
    level_array_(level_array),
    level_infer_(level_infer),
    m3_level_(),
    m3_infer_(),
    field_group_(field_group),
    num_fields_(cello::field_groups()->size(field_group_)),
    is_sync_child_(-1),
    is_sync_parent_(-1),
    is_mask_(-1),
    is_count_(-1),
    index_criteria_(index_criteria),
    num_criteria_(num_criteria)
{
  for (int i=0; i<3; i++) {
    m3_level_[i] = m3_level[i];
    m3_infer_[i] = m3_infer[i];
  }

  // Verify level_array and m3_level match
  // (note: using cello::config() rather than cello::hierarchy()
  // since the latter is not guaranteed to be initialized yet)
  int nd3[3] = {
    cello::config()->mesh_root_blocks[0],
    cello::config()->mesh_root_blocks[1],
    cello::config()->mesh_root_blocks[2] };
  int md3[3] = {
    cello::config()->mesh_root_size[0],
    cello::config()->mesh_root_size[1],
    cello::config()->mesh_root_size[2] };

  for (int axis=0; axis<cello::rank(); axis++) {
    ASSERT4("EnzoMethodInference()",
            "level-array axis %d level %d nd3 %d mismatch with level array size %d",
            axis,level_array,nd3[axis],m3_level_[axis],
            (pow(2,level_array)*nd3[axis]==m3_level_[axis]));
  }

  for (int axis=0; axis<cello::rank(); axis++) {
    int nb = md3[axis]/nd3[axis]; // compute block size along axis
    ASSERT5("EnzoMethodInference()",
            "2**(%d - %d)*(%d/%d) "
            "mismatch with inference array size %d",
            level_infer,level_array,md3[axis],nd3[axis],m3_infer_[axis],
            (pow(2,(level_infer-level_array))*nb == m3_infer_[axis]));
  }

  // Define field_group fields
  for (int i_f = 0; i_f<num_fields_; i_f++) {
    const std::string field_name =
      cello::field_groups()->item(field_group_,i_f);
    cello::define_field (field_name);
  }

  // Initialize default Refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  // BUG WORKAROUND: add a field to bypass hang bug
  refresh->add_field("density");

  // Create Level Array & broadcast proxy to other Simulation objects
  if (CkMyPe() == 0) {

    TRACE_ARRAY_CALL("ckNew()");
    proxy_level_array = CProxy_EnzoLevelArray::ckNew();

    TRACE_ENTRY_CALL("EnzoSimulation.p_set_level_array()");
    proxy_enzo_simulation.p_set_level_array(proxy_level_array);

  }

  // Create and initialize Block synchronization counters
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  ScalarDescr * scalar_descr_int  = cello::scalar_descr_int();

  // Initialize Block scalars

  //    child-to-parent synch counter
  is_sync_child_ = scalar_descr_sync->new_value
    ("method_inference_sync_child");

  //    parent-to-child synch counter
  is_sync_parent_ = scalar_descr_sync->new_value
    ("method_inference_sync_parent");

  //    mask defining overlapping inference arrays to create
  is_mask_ = scalar_descr_void->new_value
    ("method_inference_mask");

  //    count of inference arrays overlapping block
  is_count_ = scalar_descr_int->new_value
    ("method_inference_count");

  // Initialize Simulation synch counters

  const int num_root_blocks = nd3[0]*nd3[1]*nd3[2];

  //    root-level-blocks to simulation[0] sync counter
  enzo::simulation()->set_sync_infer_count(num_root_blocks);

  //    level array chare array element to simulation[0] sync counter
  enzo::simulation()->set_sync_infer_create(0);
}

//----------------------------------------------------------------------

void EnzoSimulation::p_set_level_array(CProxy_EnzoLevelArray proxy)
{
  TRACE_ENTRY_START("EnzoSimulation.p_set_level_array()");
  TRACE_ARRAY("setting proxy_level_array");
  proxy_level_array = proxy;
}

//----------------------------------------------------------------------

void EnzoMethodInference::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | level_base_;
  p | level_array_;
  p | level_infer_;
  PUParray(p,m3_level_,3);
  PUParray(p,m3_infer_,3);
  p | field_group_;
  p | num_fields_;
  p | is_sync_child_;
  p | is_sync_parent_;
  p | is_mask_;
  p | index_criteria_;
  p | num_criteria_;
}

//----------------------------------------------------------------------

void EnzoMethodInference::compute ( Block * block) throw()
{
  delete [] *scalar_mask_(block);
  *scalar_mask_(block) = nullptr;
  scalar_count_(block) = 0;

  TRACE_INFER_BLOCK(block,"compute()");
  const int level = block->level();

  // Initalize block scalars: sync counters and mask array
  const int nc = cello::num_children();
  sync_child_ (block).set_stop(nc+1); // +1 for self
  TRACE_SYNC(sync_child_(block),"sync_child set_stop()");
  sync_parent_(block).set_stop(1+1);  // +1 for self
  TRACE_SYNC(sync_parent_(block),"sync_parent set_stop()");
  *scalar_mask_(block) = nullptr;

  int count = 0;
  if (block->is_leaf()) {

    //    EnzoMsgInferCreate * msg = new EnzoMsgInferCreate;

    // Apply inference array creation criteria
    apply_criteria_(block);

    forward_create_array_(block,count);

  } else if (level >= 0) {

    // Call count_arrays() and merge_masks() on own Block to
    // ensure data dependency (+1 for self)
    if (0 <= level && level < level_base_) {

      count_arrays(block,0);

    } else {

      int n, ic3[3] = {0}; // not accessed
      char * mask = *scalar_mask_(block);
      merge_masks (block, n=0, mask=nullptr, ic3);
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodInference::apply_criteria_ (Block * block)
{
  TRACE_INFER_BLOCK(block,"apply_criteria()");
  // Apply criteria whether to create the inference array under
  // any cells. Return mask of inference arrays to create
  if (block->level() >= level_base_) {
    // allocate mask
    int mx,my,mz;
    std::tie(mx,my,mz) = mask_dims_(block->level());
    int mask_size = mask_allocate_(block,mx,my,mz);
    // initialize mask to 0
    char * mask = *scalar_mask_(block);
    std::fill_n(mask,mask_size,0);
    if (block->level() >= level_base_ + 1) {
      // Set mask according to local overdensity
      //      compute_overdensity_(block,mask,mx,my,mz,3.0);
      std::fill_n(mask,mask_size,1);
    }
    print_mask_(block);
  }
}

//----------------------------------------------------------------------

int EnzoMethodInference::mask_allocate_(Block * block, int mx, int my, int mz)
{
  const int m = mx*my*mz;
  char ** mask = scalar_mask_(block);
  // allocate level array mask if not already allocated
  if (*mask == nullptr) {
    *mask = new char[m];
    std::fill_n(*mask,m,0);
  }
  return m;
}

//----------------------------------------------------------------------

std::tuple<int,int,int> EnzoMethodInference::mask_dims_(int level) const
{
  int mx,my,mz;
  const int nc = cello::num_children();
  int mask_size = int(pow(nc,(level_array_ - level)));
  const int rank = cello::rank();
  int len = pow(2,(level_array_ - level));
  mx = (rank >= 1) ? len : 1;
  my = (rank >= 2) ? len : 1;
  mz = (rank >= 3) ? len : 1;
  mask_size = std::max(mask_size,1);
  ASSERT4 ("EnzoMethodInference::mask_dims_",
           "Unexpected mask dims %d %d %d != size %d",
           mx,my,mz,mask_size,
           mx*my*mz == mask_size);
  return {mx,my,mz};
}

//----------------------------------------------------------------------

void EnzoMethodInference::compute_overdensity_
(Block * block, char * mask, int nx, int ny, int nz, float threshold)
{
  Field field = block->data()->field();
  const int index_field = field.field_id("density");

  enzo_float * de = (enzo_float *) field.values(index_field);
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (index_field,&mx,&my,&mz);
  field.ghost_depth(index_field,&gx,&gy,&gz);

  // Compute average density de_avg
  int count = 0;
  enzo_float de_sum = 0.0;
  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
        int i=ix+mx*(iy + my*iz);
        ++count;
        de_sum += de[i];
      }
    }
  }
  enzo_float de_avg = de_sum/count;

  // Set mask according to overdensity
  for (int iz=gz; iz<mz-gz; iz++) {
    int kz=((iz-gz)*nz)/(mz-2*gz);
    for (int iy=gy; iy<my-gy; iy++) {
      int ky=((iy-gy)*ny)/(my-2*gy);
      for (int ix=gx; ix<mx-gx; ix++) {
        int kx=((ix-gx)*nx)/(mx-2*gx);
        int i=ix+mx*(iy + my*iz);
        if (de[i] > threshold*de_avg) {
          int k=kx + nx*(ky+ny*kz);
          mask[k] = 1;
        }
      }
    }
  }
}

//----------------------------------------------------------------------

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
    TRACE_ENTRY_CALL("EnzoSimulation.p_infer_set_array_count()");
    proxy_enzo_simulation[0].p_infer_set_array_count(count);

  } else {

    Index index_parent = block->index().index_parent();

    auto parent_block = enzo::block_array()[index_parent];

    if (level_base_ < level) {

      // if level > level_base, accumulate level array masks
      int ic3[3] = {0,0,0};
      block->index().child(level,ic3,ic3+1,ic3+2);
      int mx,my,mz;
      std::tie(mx,my,mz) = mask_dims_(block->level());
      const int n = mx*my*mz;
      char ** mask = scalar_mask_(block);

      // merge mask into parent's
      TRACE_ENTRY_CALL("EnzoBlock.p_method_infer_merge_masks()");
      parent_block.p_method_infer_merge_masks (n,*mask,ic3);

    } else { // 0 < level <= level_base

      // count number of inference arrays if this is the base level
      if (level == level_base_) {
        count = count_arrays_(block);
      }

      // contribute inference array count to parent's count
      TRACE_ENTRY_CALL("EnzoBlock.p_method_infer_count_arrays()");
      parent_block.p_method_infer_count_arrays (count);
    }
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_merge_masks (int n, char *mask, int ic3[3])
{
  TRACE_ENTRY_START("EnzoBlock.p_method_infer_merge_masks()");
  // Entry method for merging mask arrays

  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

  // Transfer control to EnzoMethodInference
  method->merge_masks(this,n,mask,ic3);
}

//----------------------------------------------------------------------

void EnzoMethodInference::merge_masks
(Block * block, int n_in, char *mask_in, int ic3[3])
{
  char buffer[80];
  sprintf(buffer,"merge_masks() %s",block->name().c_str());
  TRACE_INFER_BLOCK (block,buffer);

  const int level = block->level();
  
  // Merge masks here only if source is child block (not +1 sync for self)
  if (n_in > 0) {

    // allocate this block's mask (if not allocated yet)
    int mx,my,mz;
    std::tie(mx,my,mz) = mask_dims_(block->level());
    mask_allocate_(block,mx,my,mz);
    char * mask = *scalar_mask_(block);

    if (level < level_array_) {

      merge_masks_2to1_ (mask, mask_in, level, ic3);

    } else { // (level >= level_array)

      merge_masks_1to1_ (mask, mask_in, level);
    }
  }

  // Check if done
  TRACE_SYNC(sync_child_(block),"sync_child next()");
  if (sync_child_(block).next()) {
    // Recurse on block forwarding to root
    forward_create_array_ (block,0);
  }
}

//----------------------------------------------------------------------

void EnzoMethodInference::merge_masks_2to1_
(char * mask, char * mask_in, int level, int ic3[3])
{
  // get my mask size
  int mx,my,mz;
  std::tie(mx,my,mz) = mask_dims_(level);

  // get child's mask size
  const int mcx = std::max(mx/2,1);
  const int mcy = std::max(my/2,1);
  const int mcz = std::max(mz/2,1);

  // copy child's mask into my mask corner
  for (int icz=0; icz<mcz; icz++) {
    int iz = ic3[2]*mcz+icz;
    for (int icy=0; icy<mcy; icy++) {
      int iy = ic3[1]*mcy+icy;
      for (int icx=0; icx<mcx; icx++) {
        int ix = ic3[0]*mcx+icx;

        int i  = ix  + mx *(iy  + my *iz);
        int ic = icx + mcx*(icy + mcy*icz);

        mask[i] = mask_in[ic];
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodInference::merge_masks_1to1_
(char * mask, char * mask_in, int level)
{
  // get my mask size
  int mx,my,mz;
  std::tie(mx,my,mz) = mask_dims_(level);

  // merge child's mask into my mask
  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
        int i  = ix  + mx *(iy  + my *iz);
        // merge mask values
        mask[i] = std::max(mask[i],mask_in[i]);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_infer_set_array_count(int count)
{
  TRACE_ENTRY_START("EnzoSimulation.p_infer_set_array_count()");

  infer_count_arrays_ += count;
  char buffer[80];
  sprintf(buffer,"EnzoSimulation::p_infer_count_arrays %d (%d)",count,infer_count_arrays_);
  TRACE_INFER(buffer);
  TRACE_SYNC(sync_infer_count_,"sync_infer_count next()");
  if (sync_infer_count_.next()) {
    sync_infer_create_.set_stop(infer_count_arrays_ + 1);
    TRACE_SYNC(sync_infer_create_,"sync_infer_create set_stop()");
    sync_infer_done_.set_stop(infer_count_arrays_);
    TRACE_SYNC(sync_infer_done_,"sync_infer_done set_stop()");
    infer_count_arrays_ = 0;
    infer_check_create_();
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_infer_array_created()
// Count number of inference arrays (level array elements) created.
{
  TRACE_ENTRY_START("EnzoSimulation.p_infer_array_created()");
  infer_check_create_();
}

//----------------------------------------------------------------------
void EnzoSimulation::infer_check_create_()
{
  TRACE_INFER_SIMULATION("create arrays CHECK");
  TRACE_SYNC(sync_infer_create_,"sync_infer_create next()");
  if (sync_infer_create_.next()) {

    TRACE_SYNC(sync_infer_create_,"sync_infer_create check 1()");
    if (sync_infer_create_.stop() == 1) {

      // No inference arrays!  Exit
      enzo::block_array().p_method_infer_exit();

    } else {
      // All level array elements have now been created
      TRACE_INFER_SIMULATION("create arrays DONE");
      TRACE_ARRAY_CALL("doneInserting()");
      proxy_level_array.doneInserting();

      // Have level array elements request data from overlapping Blocks
      TRACE_ENTRY_CALL("EnzoLevelArray.p_request_data()");
      TRACE_ARRAY_CALL("p_request_data()");

      proxy_level_array.p_request_data();
    }
    sync_infer_create_.reset();
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::p_request_data()
{
  TRACE_ENTRY_START("EnzoLevelArray.p_request_data()");
  TRACE_ARRAY_START("p_request_data()");
  char buffer[80];
  int i3[3];
  thisIndex.values(i3);
  sprintf(buffer,"EnzoMethodInference::request_data %d %d %d %s",
          i3[0],i3[1],i3[2],field_group_.c_str());
  TRACE_INFER (buffer);

  // Get root blocking size
  int ndx,ndy,ndz;
  cello::hierarchy()->root_blocks(&ndx,&ndy,&ndz);

  // Get AMR array (vs tree) index from level array index
  int il3[3] = { thisIndex[0], thisIndex[1], thisIndex[2] };
  // shift to level_root = 0
  int iax = il3[0] >> level_array_;
  int iay = il3[1] >> level_array_;
  int iaz = il3[2] >> level_array_;

  // Determine containing base-level block

  Index index;
  index.set_array(iax,iay,iaz);
  for (int level=1; level<=level_base_; level++) {
    int tree_shift = (level_array_ - level);
    int ic3[3] = {
      (il3[0] >> tree_shift) & 1,
      (il3[1] >> tree_shift) & 1,
      (il3[2] >> tree_shift) & 1};
    index.push_child(ic3[0],ic3[1],ic3[2]);
  }

  TRACE_ENTRY_CALL("EnzoBlock.p_method_infer_request_data()");
  enzo::block_array()[index].p_method_infer_request_data(il3);
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_request_data (int ia3[3])
{
  TRACE_ENTRY_START("EnzoBlock.p_method_infer_request_data()");
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());
  method->request_data(this,ia3);
}

void EnzoMethodInference::request_data (Block * block, int ia3[3])
{
  if (block->is_leaf()) {

    // Debug output
    {
      char buffer[80];
      sprintf (buffer,"request_data leaf %d %d %d",ia3[0],ia3[1],ia3[2]);
      TRACE_INFER_BLOCK(block,buffer);
    }

    // Get and check number of fields

    Field field = block->data()->field();

    ASSERT1 ("EnzoMethodInference::request_data()",
             "Field group %s is empty or does not exist",
             field_group_.c_str(),
             (num_fields_ > 0));

    // Declare and initialize vector of vector values to send to level array

    // Compute field portion offsets and sizes for /all/ fields
    // Stored to avoid multiple calls to get_block_portion_()

    const int n_f = num_fields_;
    // field portion offsets
    std::vector<int> ofx(n_f);
    std::vector<int> ofy(n_f);
    std::vector<int> ofz(n_f);
    // field portion sizes
    std::vector<int> nfx(n_f);
    std::vector<int> nfy(n_f);
    std::vector<int> nfz(n_f);

    // compute buffer size nb

    int nb = 0;
    for (int i_f=0; i_f<n_f; i_f++) {

      const std::string field_name =
        cello::field_groups() -> item(field_group_,i_f);

      const int index_field = field.field_id (field_name);

      std::tie(ofx[i_f],ofy[i_f],ofz[i_f],nfx[i_f],nfy[i_f],nfz[i_f]) =
        get_block_portion_(block->index(), index_field, ia3);

      // Reserve storage for field value portion, plus nf[xyz] ( + 3 )

      nb += 3 + nfx[i_f]*nfy[i_f]*nfz[i_f];
    }

    // allocate buffer
    std::vector< enzo_float> buffer_values;
    buffer_values.resize(nb);

    // copy field portions

    if (block->level() > level_infer_) {

      // restrict fields before sending
      // (only if resolution higher than inference array's)
      INCOMPLETE("restrict");

    } else {

      // copy fields data

      int i_b = 0; // buffer index

      for (int i_f=0; i_f<n_f; i_f++) {

        const std::string field_name =
          cello::field_groups() -> item(field_group_,i_f);

        enzo_float * field_values = (enzo_float *)field.values(field_name);

        const int index_field = field.field_id (field_name);

        int mx,my,mz;
        field.dimensions (index_field,&mx,&my,&mz);

        // First copy of[xyz] and nf[xyz] indices to buffer
        // (to avoid having to recompute)

        buffer_values[i_b++] = nfx[i_f];
        buffer_values[i_b++] = nfy[i_f];
        buffer_values[i_b++] = nfz[i_f];

        // Then copy field portion values to buffer
        for (int iz=0; iz<nfz[i_f]; iz++) {
          const int ifz = ofz[i_f] + iz;
          for (int iy=0; iy<nfy[i_f]; iy++) {
            const int ify = ofy[i_f] + iy;
            for (int ix=0; ix<nfx[i_f]; ix++) {
              const int ifx = ofx[i_f] + ix;
              const int i_f = ifx + mx*(ify+my*ifz);
              buffer_values[i_b++] = field_values[i_f];
            }
          }
        }
      } // for i_f
      ASSERT2("EnzoMethodInference::request_data",
              "Mismatch betwen expected %d and actual %d buffer size",
              i_b,nb,
              (i_b == nb));
    }

    Index3 index3(ia3[0],ia3[1],ia3[2]);

    TRACE_ENTRY_CALL("EnzoLevelArray.p_transfer_data()");
    TRACE_INFER_BLOCK(block,"Calling p_transfer_data()");
    TRACE_ARRAY_CALL("p_transfer_data()");
    proxy_level_array[index3].p_transfer_data
      (block->index(),nb,buffer_values.data());

  } else { // not leaf

    // forward request to all intersecting child blocks
    ItChild it_child(cello::rank());
    int ic3[3];
    while (it_child.next(ic3)) {
      Index index_child = block->index().index_child(ic3);
      if (block_intersects_array_(index_child,ia3)) {
        TRACE_ENTRY_CALL("EnzoBlock.p_method_infer_request_data()");
        enzo::block_array()[index_child].p_method_infer_request_data(ia3);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::p_transfer_data
(Index index, int nb, enzo_float * buffer_values )
{
  TRACE_ARRAY_START("p_transfer_data()");
  TRACE_ENTRY_START("EnzoLevelArray.p_transfer_data()");
  // // Loop through fields
  int i_b = 0;
  const int level = index.level();
  
  for (int i_f=0; i_f<num_fields_; i_f++) {

    const int nbx = buffer_values[i_b++];
    const int nby = buffer_values[i_b++];
    const int nbz = buffer_values[i_b++];

    const std::string field_name =
      cello::field_groups() -> item(field_group_,i_f);

    const int index_field = cello::field_descr()->field_id (field_name);

    // Get initial field portion and final array values
    enzo_float * field = buffer_values + i_b;
    enzo_float * array = field_values_[i_f].data();

    const int na = nax_*nay_*naz_;
    enzo_float * a_c = field;
    enzo_float * a_f = nullptr;
    const int rank = cello::rank();
    const int rx = (rank >= 1) ? 2 : 1;
    const int ry = (rank >= 2) ? 2 : 1;
    const int rz = (rank >= 3) ? 2 : 1;
    int ecx = (rank >= 1) ? 1 : 0;
    int ecy = (rank >= 2) ? 1 : 0;
    int ecz = (rank >= 3) ? 1 : 0;
    int efx = (rank >= 1) ? 1 : 0;
    int efy = (rank >= 2) ? 1 : 0;
    int efz = (rank >= 3) ? 1 : 0;

    if (level < level_infer_) {
      // Interpolate
      // initialize coarse size
      int ncx,ncy,ncz;
      int nfx,nfy,nfz;
      for (int l=level; l<level_infer_; l++) {
        const bool is_first = (l == level);
        const bool is_last  = (l == level_infer_ - 1);
        ncx = is_first ? nbx : (ncx-2*ecx)*rx+2*ecx;
        ncy = is_first ? nby : (ncy-2*ecy)*ry+2*ecy;
        ncz = is_first ? nbz : (ncz-2*ecz)*rz+2*ecz;
        a_c = is_first ? field : a_f;

        nfx = is_last ? (ncx-2*efx)*rx : (ncx-2*efx)*rx+2*efx;
        nfy = is_last ? (ncy-2*efy)*ry : (ncy-2*efy)*ry+2*efy;
        nfz = is_last ? (ncz-2*efz)*rz : (ncz-2*efz)*rz+2*efz;
        a_f = is_last ? array : new enzo_float [nfx*nfy*nfz];
        efx = is_last ? 0 : efx;
        efy = is_last ? 0 : efy;
        efz = is_last ? 0 : efz;
        interpolate_(a_f,nfx,nfy,nfz,nfx,nfy,nfz,efx,efy,efz,
                     a_c,ncx,ncy,ncz,ncx,ncy,ncz,ecx,ecy,ecz);
        if (!is_last) delete [] a_f;
      }
    } else if (level > level_infer_) {
      // Restrict
      for (int l=level_infer_; l<level; l++) {
        INCOMPLETE2 ("p_transfer_data() restrict %d to %d",l,l-1);
      }
    } else {
      // Copy
      int l = level;
      int mdx = nbx+2*efx;
      int mdy = nby+2*efy;
      int mdz = nbz+2*efz;
      copy_(a_f,nbx,nby,nbz,nbx,nby,nbz,0,0,0,
            a_c,mdx,mdy,mdz,nbx,nby,nbz,efx,efy,efz);
    }

    // advance pointer to next field
    i_b += nbx*nby*nbz;

  }

  ASSERT2("EnzoLevelArray::p_transfer_data()",
          "Mismatch betwen expected %d and actual %d buffer size",
          i_b,nb,
          (i_b == nb));

  // Check if this is the last block portion expected by tracking volumes

  int level_volume = std::max(0,index.level() - level_array_);
  float volume_portion = pow(1.0/cello::num_children(),level_volume);
  volume_ratio_ += volume_portion;

  ASSERT1 ("EnzoLevelArray::p_transfer_data()",
           "volume_ratio_ %g not between 0.0 and 1.0 as expected",
           volume_ratio_,
           (0.0 <= volume_ratio_) && (volume_ratio_ <= 1.0));

  if (volume_ratio_ == 1.0) {
    volume_ratio_ = 0.0;
    // When done, call inference (note safe to compare float with constant
    // since guaranteed no roundoff error)
    apply_inference();
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::apply_inference()
{
  // Apply inference method

  // ...

  // Synchronize when done
  proxy_enzo_simulation[0].p_infer_done();
}

//----------------------------------------------------------------------

void EnzoSimulation::p_infer_done()
{
  TRACE_SYNC(sync_infer_done_,"sync_infer_done next()");
  if (sync_infer_done_.next()) {
    sync_infer_done_.reset();
    TRACE_SYNC(sync_infer_done_,"sync_infer_done reset()");

    TRACE_INFER_SIMULATION("inference DONE");

    TRACE_ARRAY_CALL("ckDestroy()");
    proxy_level_array.ckDestroy();

    // Have EnzoBlock elements exit
    enzo::block_array().p_method_infer_exit();
  }
}

//----------------------------------------------------------------------

std::tuple<int,int,int,int,int,int>
EnzoMethodInference::get_block_portion_
(Index index, int index_field, int ia3[3])
{
  int gx,gy,gz;
  cello::field_descr()->ghost_depth(index_field,&gx,&gy,&gz);

  const int rank = cello::rank();

  // Get block size (nx,ny,nz)
  int nx,ny,nz;

  {
    int ndx,ndy,ndz;
    int mdx,mdy,mdz;
    cello::hierarchy()->root_blocks(&ndx,&ndy,&ndz);
    cello::hierarchy()->root_size  (&mdx,&mdy,&mdz);
    nx = mdx/ndx;
    ny = mdy/ndy;
    nz = mdz/ndz;
  }

  // Determine reduced size if restricting
  const int num_restrict   = std::max(index.level() - level_infer_,0);
  const int ratio_restrict = pow(2,num_restrict);

  // Determine reduced size if multiple inference arrays per block
  const int num_infer   = std::max(level_array_ - index.level(),0);
  const int ratio_infer = pow(2,num_infer);

  // get offset into inference array
  int kbx=0;
  int kby=0;
  int kbz=0;
  if (ratio_infer > 1) {
    int mask = (ratio_infer - 1);
    kbx = (ia3[0] & mask);
    kby = (ia3[1] & mask);
    kbz = (ia3[2] & mask);
  }

  // get extra padding amount
  const int ex = (rank >= 1)? 1 : 0;
  const int ey = (rank >= 2)? 1 : 0;
  const int ez = (rank >= 3)? 1 : 0;

  // compute field portion size
  int nbx,nby,nbz;

  const int ratio = ratio_restrict*ratio_infer;

  nbx = (rank >= 1) ? nx / ratio + 2*ex : 1;
  nby = (rank >= 2) ? ny / ratio + 2*ey : 1;
  nbz = (rank >= 3) ? nz / ratio + 2*ez : 1;

  // compute field portion offset
  int obx,oby,obz;

  obx = kbx*nx/ratio + gx - ex;
  oby = kby*ny/ratio + gy - ey;
  obz = kbz*nz/ratio + gz - ez;

  return {obx,oby,obz,nbx,nby,nbz};
}

//----------------------------------------------------------------------

void EnzoLevelArray::coarsen_
(enzo_float * a_c,
 int mcx, int mcy, int mcz, int ncx, int ncy, int ncz, int ecx, int ecy, int ecz,
 const enzo_float * a_f,
 int mfx, int mfy, int mfz, int nfx, int nfy, int nfz, int efx, int efy, int efz)
{
  const int dx = 1;
  const int dy = mcx;
  const int dz = mcx*mcz;
  const int d000 = 0;
  const int d001 = dx;
  const int d010 = dy;
  const int d011 = dx+dy;
  const int d100 = dz;
  const int d101 = dx+dz;
  const int d110 = dx+dy;
  const int d111 = dx+dy+dz;
  if (cello::rank() == 1) {
  } else if (cello::rank() == 2) {
  } else if (cello::rank() == 3) {
    for (int icz=0; icz<mcz; icz++) {
      for (int icy=0; icy<mcy; icy++) {
        for (int icx=0; icx<mcx; icx++) {
          int i_f=icx*2 + mfx*(icy*2 + mfy*icz*2);
          int i_c=icx   + mcx*(icy   + mcy*icz);
          a_c[i_c] = 0.125*(a_f[i_f+d000] + a_f[i_f+d001] + a_f[i_f+d010] + a_f[i_f+d011] +
                            a_f[i_f+d100] + a_f[i_f+d101] + a_f[i_f+d110] + a_f[i_f+d111]);
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::interpolate_
(enzo_float * a_f,
 int mfx, int mfy, int mfz, int nfx,int nfy, int nfz, int efx, int efy, int efz,
 const enzo_float * a_c,
 int mcx, int mcy, int mcz, int ncx,int ncy, int ncz, int ecx, int ecy, int ecz)
{
  const int dx = 1;
  const int dy = mcx;
  const int dz = mcx*mcz;
  const int d000 = 0;
  const int d001 = dx;
  const int d010 = dy;
  const int d011 = dx+dy;
  const int d100 = dz;
  const int d101 = dx+dz;
  const int d110 = dx+dy;
  const int d111 = dx+dy+dz;

  const enzo_float w[2] = { 0.75, 0.25 };

  if (cello::rank() == 1) {
    for (auto ifx=0; ifx<nfx; ifx++) {
      auto iex = ifx+(1-efx);
      auto icx = iex/2;
      auto wx0 = w[iex&1];
      auto wx1 = 1.0-wx0;
      auto i_f = ifx;
      auto i_c = icx;
      a_f[i_f] = 0.5*(wx0*a_c[i_c+d000] +
                      wx1*a_c[i_c+d001]);
    }
  } else if (cello::rank() == 2) {
    for (auto ify=0; ify<nfy; ify++) {
      auto iey = ify+(1-efy);
      auto icy = iey/2;
      auto wy0 = w[iey&1];
      auto wy1 = 1.0-wy0;
      for (auto ifx=0; ifx<nfx; ifx++) {
        auto iex = ifx+(1-efx);
        auto icx = iex/2;
        auto wx0 = w[iex&1];
        auto wx1 = 1.0-wx0;
        auto i_f = ifx + mfx*ify;
        auto i_c = icx + mcx*icy;
        a_f[i_f] = 0.25*(wy0*wx0*a_c[i_c+d000] +
                         wy0*wx1*a_c[i_c+d001] +
                         wy1*wx0*a_c[i_c+d010] +
                         wy1*wx1*a_c[i_c+d011]);
      }
    }
  } else if (cello::rank() == 3) {
    for (auto ifz=0; ifz<nfz; ifz++) {
      auto iez = ifz+(1-efz);
      auto icz = iez/2;
      auto wz0 = w[iez&1];
      auto wz1 = 1.0-wz0;
      for (auto ify=0; ify<nfy; ify++) {
        auto iey = ify+(1-efy);
        auto icy = iey/2;
        auto wy0 = w[iey&1];
        auto wy1 = 1.0-wy0;
        for (auto ifx=0; ifx<nfx; ifx++) {
          auto iex = ifx+(1-efx);
          auto icx = iex/2;
          auto wx0 = w[iex&1];
          auto wx1 = 1.0-wx0;
          auto i_f = ifx + mfx*(ify + mfy*ifz);
          auto i_c = icx + mcx*(icy + mcy*icz);
          a_f[i_f] = 0.125*(wz0*wy0*wx0*a_c[i_c+d000] +
                            wz0*wy0*wx1*a_c[i_c+d001] +
                            wz0*wy1*wx0*a_c[i_c+d010] +
                            wz0*wy1*wx1*a_c[i_c+d011] +
                            wz1*wy0*wx0*a_c[i_c+d100] +
                            wz1*wy0*wx1*a_c[i_c+d101] +
                            wz1*wy1*wx0*a_c[i_c+d110] +
                            wz1*wy1*wx1*a_c[i_c+d111]);
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::copy_
(enzo_float * a_dst,
 int mdx, int mdy, int mdz, int ndx,int ndy, int ndz, int edx, int edy, int edz,
 const enzo_float * a_src,
 int msx, int msy, int msz, int nsx,int nsy, int nsz, int esx, int esy, int esz)
{
  if (cello::rank() == 1) {
    for (int idx=0; idx<ndx; idx++) {
      const int isx = idx;
      const int id= idx;
      const int is= isx;
      a_dst[id] = a_src[is];
    }
  } else if (cello::rank() == 2) {
    for (int idy=0; idy<ndy; idy++) {
      const int isy = idy;
      for (int idx=0; idx<ndx; idx++) {
        const int isx = idx;
        const int id= idx+mdx*idy;
        const int is= isx+msx*isy;
        a_dst[id] = a_src[is];
      }
    }
  } else if (cello::rank() == 3) {
    for (int idz=0; idz<ndz; idz++) {
      const int isz = idz;
      for (int idy=0; idy<ndy; idy++) {
        const int isy = idy;
        for (int idx=0; idx<ndx; idx++) {
          const int isx = idx;
          const int id= idx+mdx*(idy+mdy*idz);
          const int is= isx+msx*(isy+msy*isz);
          a_dst[id] = a_src[is];
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_count_arrays (int count)
{
  TRACE_ENTRY_START("EnzoBlock.p_method_infer_count_arrays()");
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

  method->count_arrays(this,count);
}

void EnzoMethodInference::count_arrays (Block * block, int count)
{
  char buffer[80];
  scalar_count_(block) += count;
  sprintf (buffer,"count_arrays() count %d (%d)\n",count,scalar_count_(block));
  TRACE_INFER_BLOCK (block,buffer);
  TRACE_SYNC(sync_child_(block),"sync_child next()");

  if (sync_child_(block).next()) {
    // Recurse on block forwarding to root
    forward_create_array_ (block,scalar_count_(block));
  }
}

//----------------------------------------------------------------------

int EnzoMethodInference::count_arrays_ (Block * block) const
{
  int mx,my,mz;
  std::tie(mx,my,mz) = mask_dims_(block->level());
  const int n=mx*my*mz;
  const char * mask = *scalar_mask_(block);
  int count = 0;
  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
        const int i=ix+mx*(iy+my*iz);
        if (mask[i]) ++count;
      }
    }
  }
  return count;
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_exit()
{
  TRACE_ENTRY_START("EnzoBlock.p_method_infer_exit()");
  TRACE_INFER_BLOCK(this,"p_method_infer_exit()");
  compute_done();
}

//======================================================================

void EnzoMethodInference::create_level_arrays_ (Block * block)
{
  TRACE_INFER_BLOCK(block,"create_level_arrays_()");
  ASSERT2 ("EnzoMethodInference::create_level_arrays_",
           "Block %s not in expected refinement level %d",
           block->name().c_str(),block->level(),
           block->level() == level_base_);
  int mx,my,mz;
  std::tie(mx,my,mz) = mask_dims_(block->level());
  int m = mx*my*mz;
  char * mask = *scalar_mask_(block);
  int nax,nay,naz;
  level_array_dims_(&nax,&nay,&naz);
  if (mask != nullptr) {
    int i3[3];
    block->index().index_level (i3, block->level());
    for (int kz=0; kz<mz; kz++) {
      for (int ky=0; ky<my; ky++) {
        for (int kx=0; kx<mx; kx++) {
          const int k=kx+mx*(ky + my*kz);
          if (mask[k]) {
            Index3 index(i3[0]*mx+kx,
                         i3[1]*my+ky,
                         i3[2]*mz+kz);
            TRACE_ARRAY_CALL("insert()");
            proxy_level_array[index].insert
              (field_group_,level_base_, level_array_,level_infer_,nax,nay,naz);
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
  cello::hierarchy()->root_blocks(mx,my,mz,level_array_);
}

bool EnzoMethodInference::block_intersects_array_(Index index, int ia3[3])
{
  // get lower-index of block in level_array level
  int ib3[3];
  index.index_level(ib3,level_array_);
  // get upper-index
  const int r = pow(2,(level_array_ - index.level()));
  int lx = (ib3[0] <= ia3[0] && ia3[0] < ib3[0]+r);
  int ly = (ib3[1] <= ia3[1] && ia3[1] < ib3[1]+r);
  int lz = (ib3[2] <= ia3[2] && ia3[2] < ib3[2]+r);
  bool intersects = lx && ly && lz;
  return intersects;
}
