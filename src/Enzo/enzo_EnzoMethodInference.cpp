// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodInference.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-07-12
/// @brief    Implements the EnzoMethodInference class

#include "cello.hpp"
#include "enzo.hpp"

#define TRACE_INFER

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
#else
#  define TRACE_INFER(MSG) /* ... */
#  define TRACE_INFER_BLOCK(B,MSG) /* ... */
#  define TRACE_INFER_SIMULATION(SIM,MSG) /* ... */
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
    is_sync_child_(-1),
    is_sync_parent_(-1),
    is_mask_(-1),
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
    ASSERT2("EnzoMethodInference()",
            "level-array level %d mismatch with level array size %d",
            level_array,m3_level_[axis],
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
  const int num_fields = cello::field_groups()->size(field_group_);
  for (int i_f = 0; i_f<num_fields; i_f++) {
    const std::string field_name = cello::field_groups()->item(field_group_,i_f);
    cello::define_field (field_name);
  }

  // Initialize default Refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  // BUG WORKAROUND: add a field to bypass hang bug
  refresh->add_field("density");

  // Create Level Array & broadcast proxy to other Simulation objects
  if (CkMyPe() == 0) {

    proxy_level_array = CProxy_EnzoLevelArray::ckNew();

    proxy_enzo_simulation.p_set_level_array(proxy_level_array);

  }

  // Create and initialize Block synchronization counters
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();

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

  // Initialize Simulation synch counters

  const int num_root_blocks = nd3[0]*nd3[1]*nd3[2];

  //    root-level-blocks to simulation[0] sync counter
  enzo::simulation()->set_sync_infer_count(num_root_blocks);

  //    level array chare array element to simulation[0] sync counter
  enzo::simulation()->set_sync_infer_create(0);
}

//----------------------------------------------------------------------

void EnzoSimulation::p_set_level_array(CProxy_EnzoLevelArray proxy)
{ proxy_level_array = proxy; }

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
  p | is_sync_child_;
  p | is_sync_parent_;
  p | is_mask_;
  p | index_criteria_;
  p | num_criteria_;
}

//----------------------------------------------------------------------

void EnzoMethodInference::compute ( Block * block) throw()
{
  TRACE_INFER_BLOCK(block,"compute()");
  const int level = block->level();

  CkPrintf ("DEBUG_INFER %d num_blocks (level_base) %d\n",CkMyPe(),
            cello::simulation()->num_blocks_level(level_base_));

  // Initalize block scalars: sync counters and mask array
  const int nc = cello::num_children();
  sync_child_ (block).set_stop(nc+1); // +1 for self
  sync_parent_(block).set_stop(1+1);  // +1 for self
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
  int mask_size = 0;
  char ** mask = scalar_mask_(block);
  if (block->level() >= level_base_) {
    // allocate mask
    mask_size = mask_allocate_(block);
    // initialize mask
    char value = (block->level() >= level_base_ + 1) ? 1 : 0;
    //int value = 1;
    std::fill_n(*mask,mask_size,value);
  }
  CkPrintf ("mask_size %d\n",mask_size);
  CkPrintf ("Block mask %p\n",*scalar_mask_(block));
}

//----------------------------------------------------------------------

int EnzoMethodInference::mask_allocate_(Block * block)
{
  int mx,my,mz;
  std::tie(mx,my,mz) = mask_dims_(block->level());
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
      parent_block.p_method_infer_merge_masks (n,*mask,ic3);

    } else { // 0 < level <= level_base

      // count number of inference arrays if this is the base level
      if (level == level_base_) {
        count = count_arrays_(block);
      }

      // contribute inference array count to parent's count
      parent_block.p_method_infer_count_arrays (count);
    }
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_merge_masks (int n, char *mask, int ic3[3])
{
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
    mask_allocate_(block);
    char * mask = *scalar_mask_(block);

    if (level < level_array_) {

      merge_masks_2to1_ (mask, mask_in, level, ic3);

    } else { // (level >= level_array)

      merge_masks_1to1_ (mask, mask_in, level);
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
  char buffer[80];
  sprintf(buffer,"p_infer_count_arrays %d",count);

  infer_count_arrays_ += count;
  if (sync_infer_count_.next()) {
    sync_infer_create_.set_stop(infer_count_arrays_ + 1);
    infer_count_arrays_ = 0;
    infer_check_create_();
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_infer_array_created()
// Count number of inference arrays (level array elements) created.
{
  TRACE_INFER_SIMULATION("p_infer_array_created()");
  sync_infer_create_.print("array_created");
  infer_check_create_();
}

//----------------------------------------------------------------------
void EnzoSimulation::infer_check_create_()
{
  if (sync_infer_create_.next()) {
    // All level array elements have now been created
    TRACE_INFER_SIMULATION("create arrays DONE");
    proxy_level_array.doneInserting();

    // Have level array elements request data from overlapping Blocks
    proxy_level_array.p_request_data();
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::p_request_data()
{
  char buffer[80];
  int i3[3];
  thisIndex.values(i3);
  sprintf(buffer,"EnzoMethodInference::request_data %d %d %d %s",
          i3[0],i3[1],i3[2],field_group_.c_str());
  TRACE_INFER (buffer);

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // Get root blocking size
  int ndx,ndy,ndz;
  cello::hierarchy()->root_blocks(&ndx,&ndy,&ndz);

  // Get AMR array (vs tree) index from level array index
  int il3[3] = { thisIndex[0], thisIndex[1], thisIndex[2] };
  // shift to level_root = 0
  int iax = il3[0] >> level_array_;
  int iay = il3[1] >> level_array_;
  int iaz = il3[2] >> level_array_;
  CkPrintf ("DEBUG_INFER p_request_data il3 %d %d %d\n",iax,iay,iaz);

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

  enzo::block_array()[index].p_method_infer_request_data(il3);
}

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_request_data (int ia3[3])
{
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());
  method->request_data(this,ia3);
}

void EnzoMethodInference::request_data (Block * block, int ia3[3])
{
  if (block->is_leaf()) {
    {
      char buffer[80];
      sprintf (buffer,"request_data leaf %d %d %d",ia3[0],ia3[1],ia3[2]);
      TRACE_INFER_BLOCK(block,buffer);
    }

    Field field = block->data()->field();

    // Get and check number of fields
    const int num_fields = cello::field_groups()->size(field_group_);
    ASSERT1 ("EnzoMethodInference::request_data()",
             "Field group %s is empty or does not exist",
             field_group_.c_str(),
             (num_fields > 0));

    // determine whether to refine (dl < 0) copy (dl==0) or restrict (dl > 0)
    const int dl = block->level() - level_infer_;
    const int nr = dl > 0 ? dl : 0; // restrict count
    const int rp = dl < 0 ? -dl : 0; // prolong count
    const int r = pow(2,nr); // restrict ratio 1, 2, 4, etc.
    const int p = pow(2,rp); // prolong ratio 1, 2, 4, etc.

    // Declare and initialize vector of vector values to send to level array

    // get block size
    int nb3[3];
    field.size(nb3,nb3+1,nb3+2);
    // compute array size (smaller if restricting first)
    const int nx = nb3[0]/r;
    const int ny = nb3[1]/r;
    const int nz = nb3[2]/r;


    // Allocate field data, including one ghost zone
    const int nb = (nx+2)*(ny+2)*(nz+2);
    std::vector< enzo_float> buffer_values;
    buffer_values.resize(num_fields*nb);

    if (block->level() > level_infer_) {
      // restrict fields before sending
      INCOMPLETE("restrict");
    } else {
      // copy fields (prolong done at receiver)
      for (int i_f=0; i_f<num_fields; i_f++) {
        const std::string field_name = cello::field_groups()->item(field_group_,i_f);
        enzo_float * field_values = (enzo_float *)field.values(field_name);
        const int index_field = field.field_id (field_name);
        int mx,my,mz;
        field.dimensions(index_field,&mx,&my,&mz);
        const int gx = (mx-nx)/2;
        const int gy = (my-ny)/2;
        const int gz = (mz-nz)/2;
        for (int iz=0; iz<nz; iz++) {
          const int ifz = iz+gz;
          for (int iy=0; iy<ny; iy++) {
            const int ify = iy+gy;
            for (int ix=0; ix<nx; ix++) {
              const int ifx = ix+gx;
              const int i_a = i_f*nb + ix + nx*(iy+ny*iz);
              const int i_f = ifx + mx*(ify+my*ifz);
              buffer_values[i_a] = field_values[i_f];
            }
          }
        }
      }
    }

    Index3 index3(ia3[0],ia3[1],ia3[2]);
    CkPrintf ("DEBUG_INFER p_transfer_data %s\n",block->name().c_str());
    proxy_level_array[index3].p_transfer_data(block->index(),num_fields*nb,buffer_values.data());
    // else if not leaf, forward request to children overlapping array
  } else {
    TRACE_INFER_BLOCK(block,"request_data not leaf");
    ItChild it_child(cello::rank());
    int ic3[3];
    Index index = block->index();
    // forward request to all intersecting child blocks
    while (it_child.next(ic3)) {
      Index index_child = index.index_child(ic3);
      if (block_intersects_array_(index_child,ia3)) {
        enzo::block_array()[index_child].p_method_infer_request_data(ia3);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::p_transfer_data (Index index, int nf, enzo_float * field_data_list )
{
  // Copy incoming field data to inference arrays, possibly with linear interpolation

  // Determine array and field sizes and counts

  const int num_fields = cello::field_groups()->size(field_group_);

  // get domain size and blocking to get block size
  int ndx,ndy,ndz;
  int mdx,mdy,mdz;
  cello::hierarchy()->root_blocks(&ndx,&ndy,&ndz);
  cello::hierarchy()->root_size  (&mdx,&mdy,&mdz);
  // block size
  const int nbx = mdx/ndx;
  const int nby = mdy/ndy;
  const int nbz = mdz/ndz;
  const int nb = nbx*nby*nbz;
  
  const int iax = thisIndex[0];
  const int iay = thisIndex[1];
  const int iaz = thisIndex[2];
  CkPrintf ("DEBUG_INFEX p_transfer_data array %d %d %d size %d fields %d block size %d %d %d\n",
            iax,iay,iaz,nf,num_fields,nbx,nby,nbz);
  ASSERT5("EnzoLevelArray::p_transfer_data()",
          "array size %d does not match field count %d times block size %d * %d * %d",
          nf,num_fields,nbx,nby,nbz,
          (nf == num_fields*(nbx+2)*(nby+2)*(nbz+2)));
  // array size
  const int na = nax_*nay_*naz_;

  // portion of block

  // ma: number of inference arrays per block
  // kaxyz: index of inferece array in block

  // compute offset 0 <= (kax,kay,kaz) < ma into field data
  bool la = level_array_ > index.level();

  const int ma = la ? (1 << (level_array_ - index.level())) : 1;
  unsigned mask_a = (ma  - 1);
  CkPrintf ("DEBUG_INFER ma %d mask %x\n",ma,mask_a);
  const int kax = iax & mask_a;
  const int kay = iay & mask_a;
  const int kaz = iaz & mask_a;
  CkPrintf ("DEBUG_INFER kaxyz %d %d %d\n",kax,kay,kaz);
  const int dax = kax*nax_/ma;
  const int day = kay*nay_/ma;
  const int daz = kaz*naz_/ma;
  CkPrintf ("DEBUG_INFER daxyz %d %d %d\n",dax,day,daz);
  // compute offset 0 <= (kfx,kfy,kfz) < mf into array
  //  const int kfx = ifx & mask_f;
  //  const int kfy = ify & mask_f;
  //  const int kfz = ifz & mask_f;

  bool lf = index.level() > level_array_;
  const int mf = lf ? (1 << (index.level() - level_array_)) : 1;
  unsigned mask_f = (mf  - 1);
  int t3[3];
  index.tree(t3,t3+1,t3+2,index.level());
  const int kfx = lf ? t3[0] & mask_f : 0;
  const int kfy = lf ? t3[1] & mask_f : 0;
  const int kfz = lf ? t3[2] & mask_f : 0;
  CkPrintf ("DEBUG_INFER tree 0 %x %x %x : %d %d %d\n",t3[0],t3[1],t3[2],kfx,kfy,kfz);

  CkPrintf ("DEBUG_INFER p_transfer_data() kaxyz %d %d %d mask %x/%d\n",kax,kay,kaz,mask_a,ma);
  //CkPrintf ("DEBUG_INFER p_transfer_data() kfxyz %d %d %d mask %x/%d\n",kfx,kfy,kfz,mask_f,mf);

  // Loop through fields
  for (int i_f=0; i_f<num_fields; i_f++) {

    // Field data
    enzo_float * field = field_data_list + i_f*nb;
    enzo_float * array = field_values_[i_f].data();

    // array   field_values[i_f][ia];
    // field field_data_list[i_f*nb]
    // interpolate if needed

    // Linear interpolate X times to reach level_infer_
    enzo_float * f = field;
    for (int level=index.level(); level<level_infer_; level++) {
      bool is_first = (level == index.level());
      bool is_last  = (level == level_infer_);
      enzo_float * a = (is_last) ? array : new enzo_float[na];
      const int max = nax_+dax;
      const int may = nay_+day;
      const int maz = naz_+daz;
      interpolate_(a,max,may,maz,f,nbx,nby,nbz);
      CkPrintf ("TRACE_INFER interpolate %d to %d\n",level,level+1);
    }
    if (level_infer_ == index.level()) {
      CkPrintf ("TRACE_INFER copy\n");
    }
    for (int level=level_infer_+1; level<index.level(); level++) {
      CkPrintf ("TRACE_INFER restrict %d to %d\n",level,level+1);
    }
    {
      char buffer[256];
      sprintf (buffer,"EnzoLevelArray::p_transfer_data %d %d %d  %d %g %g level %d ratio %g",
               thisIndex[0],thisIndex[1],thisIndex[2],
               nf,field_data_list[0],field_data_list[nf-1],index.level(),volume_ratio_);
      TRACE_INFER(buffer);
    }
  }

  // If this is the last block data expected,
  // compute ratio of incoming field data to update volume_ratio_
  int dl = std::max(0,index.level() - level_array_);
  volume_ratio_ += pow(1.0/cello::num_children(),dl);
  ASSERT1 ("EnzoLevelArray::p_transfer_data()",
           "volume_ratio_ %g not between 0.0 and 1.0 as expected",
           volume_ratio_,
           (0.0 <= volume_ratio_) && (volume_ratio_ <= 1.0));
  if (volume_ratio_ == 1.0) {
    CkPrintf ("TRACE_INFER %d %d %d p_transfer_data done()\n",thisIndex[0],thisIndex[1],thisIndex[2]);
    // done
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::coarsen_
(enzo_float * a_c, int mcx, int mcy, int mcz,
 enzo_float * a_f, int mfx, int mfy, int mfz)
{
  TRACE_INFER("coarsen\n");
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
(enzo_float * a_c, int mcx, int mcy, int mcz,
 enzo_float * a_f, int mfx, int mfy, int mfz)
{
  CkPrintf ("TRACE_INFER interpolate mcxyz %d %d %d  mfxyz %d %d %d\n",
              mcx,mcy,mcz,mfx,mfy,mfz);
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
    for (int icz=0; icz<mcz-1; icz++) {
      for (int icy=0; icy<mcy-1; icy++) {
        for (int icx=0; icx<mcx-1; icx++) {
          int i_f=icx*2 + mfx*(icy*2 + mfy*icz*2);
          int i_c=icx   + mcx*(icy   + mcy*icz);
          a_f[i_f] = 0.125*(a_c[i_c+d000] + a_c[i_c+d001] + a_c[i_c+d010] + a_c[i_c+d011] +
                            a_c[i_c+d100] + a_c[i_c+d101] + a_c[i_c+d110] + a_c[i_c+d111]);
        }
      }
    }
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
    // Recurse on block forwarding to root
    forward_create_array_ (block,count);
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
  int nbx,nby,nbz;
  cello::hierarchy()->root_blocks(&nbx,&nby,&nbz,level_base_);
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
