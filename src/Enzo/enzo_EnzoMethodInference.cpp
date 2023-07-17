// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodInference.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-07-12
/// @brief    Implements the EnzoMethodInference class

#include "cello.hpp"
#include "enzo.hpp"

// Used to write inference array bounds and "bubble" coordinates for plotting
// #define TRACE_INFER
// #define DEBUG_INFER

// # define PRINT_FIELD_VALUES

//----------------------------------------------------------------------

EnzoMethodInference::EnzoMethodInference
(int level_base,
 int level_array,
 int level_infer,
 const int m3_level[3],
 const int m3_infer[3],
 std::string field_group,
 float overdensity_threshold)
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
    overdensity_threshold_(overdensity_threshold)
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

    proxy_level_array = CProxy_EnzoLevelArray::ckNew();

    proxy_enzo_simulation.p_set_level_array(proxy_level_array);

  }

  // Create and initialize Block synchronization counters
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  ScalarDescr * scalar_descr_int  = cello::scalar_descr_int();

  // Initialize Block scalars

  //    child-to-parent sync counter
  is_sync_child_ = scalar_descr_sync->new_value
    ("method_inference_sync_child");

  //    parent-to-child sync counter
  is_sync_parent_ = scalar_descr_sync->new_value
    ("method_inference_sync_parent");

  //    mask defining overlapping inference arrays to create
  is_mask_ = scalar_descr_void->new_value
    ("method_inference_mask");

  //    count of inference arrays overlapping block
  is_count_ = scalar_descr_int->new_value
    ("method_inference_count");

  // Initialize Simulation sync counters

  const int num_root_blocks = nd3[0]*nd3[1]*nd3[2];

  //    root-level-blocks to simulation[0] sync counter
  enzo::simulation()->set_sync_infer_count(num_root_blocks);

  //    level array chare array element to simulation[0] sync counter
  enzo::simulation()->set_sync_infer_create(0);

}

//----------------------------------------------------------------------

void EnzoSimulation::p_set_level_array(CProxy_EnzoLevelArray proxy)
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
  p | level_array_;
  p | level_infer_;
  PUParray(p,m3_level_,3);
  PUParray(p,m3_infer_,3);
  p | field_group_;
  p | num_fields_;
  p | is_sync_child_;
  p | is_sync_parent_;
  p | is_mask_;
  p | is_count_;
  p | overdensity_threshold_;
}

//----------------------------------------------------------------------

void EnzoMethodInference::compute ( Block * block) throw()
{
  // Initalize block scalars: sync counters and mask array

  delete [] *scalar_mask_(block);
  *scalar_mask_(block) = nullptr;
  scalar_count_(block) = 0;

  const int nc = cello::num_children();
  sync_child_ (block).set_stop(nc+1); // +1 for self
  sync_parent_(block).set_stop(1+1);  // +1 for self

  const int level = block->level();

  // set random seed for models with stochastic elements
  srand( time(NULL)*(CkMyPe()+1) ); 

  if (block->is_leaf()) {

    // Apply inference array creation criteria
    apply_criteria_(block);

    forward_create_array_(block,0);

  } else if (level >= 0) {

    // Call count_arrays() or merge_masks() on this Block to
    // ensure data dependency (+1 for self)

    if (level < level_base_) {

      // call count_arrays() for levels below level_base
      count_arrays(block,0);

    } else { // level >= level_base_

      // call merge masks() for levels above and including level_base
      int n, ic3[3] = {0}; // not accessed
      char * mask = *scalar_mask_(block);
      merge_masks (block, n=0, mask=nullptr, ic3);
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodInference::apply_criteria_ (Block * block)
{
  // Apply criteria whether to create the inference array under
  // any cells. Return mask array of inference arrays to create
  if (block->level() >= level_base_) {
    // allocate mask
    int mx,my,mz;
    std::tie(mx,my,mz) = mask_dims_(block->level());
    int mask_size = mx*my*mz;
    mask_allocate_(block,mx,my,mz);
    // initialize mask to 0
    char * mask = *scalar_mask_(block);
    std::fill_n(mask,mask_size,0);
    if (block->level() >= level_array_ + 1) {
      // Set mask according to local overdensity
      compute_overdensity_(block,mask,mx,my,mz);

      //TODO: Add parameter to toggle inference mesh creation criteria
      // Set mask by checking if mean overdensity > overdensity_threshold
      //compute_mean_overdensity_(block,mask,mx,my,mz);
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodInference::mask_allocate_(Block * block, int mx, int my, int mz)
{
  // allocate level array mask if not already allocated
  char ** mask = scalar_mask_(block);
  if (*mask == nullptr) {
    const int m = mx*my*mz;
    *mask = new char[m];
    std::fill_n(*mask,m,0);
  }
}

//----------------------------------------------------------------------

std::tuple<int,int,int> EnzoMethodInference::mask_dims_(int level) const
{
  int nx,ny,nz;
  const int len = std::max(1,int(std::pow(2,(level_array_ - level))));
  const int rank = cello::rank();
  nx = (rank >= 1) ? len : 1;
  ny = (rank >= 2) ? len : 1;
  nz = (rank >= 3) ? len : 1;

  return {nx,ny,nz};
}

//----------------------------------------------------------------------

void EnzoMethodInference::compute_overdensity_
(Block * block, char * mask, int nx, int ny, int nz)
{
  Field field = block->data()->field();
  const int index_field = field.field_id("density");

  enzo_float * de = (enzo_float *) field.values(index_field);
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (index_field,&mx,&my,&mz);
  field.ghost_depth(index_field,&gx,&gy,&gz);

  // Compute average block-local density de_avg
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
  const enzo_float de_avg = de_sum/count;

  // Set mask = 1 where baryon density greater than
  // threshold times local average
  for (int iz=gz; iz<mz-gz; iz++) {
    int kz=((iz-gz)*nz)/(mz-2*gz);
    for (int iy=gy; iy<my-gy; iy++) {
      int ky=((iy-gy)*ny)/(my-2*gy);
      for (int ix=gx; ix<mx-gx; ix++) {
        int kx=((ix-gx)*nx)/(mx-2*gx);
        int i=ix + mx*(iy + my*iz);
        if (de[i] > overdensity_threshold_*de_avg) {
          int k=kx + nx*(ky + ny*kz);
          mask[k] = 1;
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodInference::compute_mean_overdensity_
(Block * block, char * mask, int nx, int ny, int nz)
{
  Field field = block->data()->field();
  const int index_field = field.field_id("density");

  enzo_float * de = (enzo_float *) field.values(index_field);
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (index_field,&mx,&my,&mz);
  field.ghost_depth(index_field,&gx,&gy,&gz);

  // Compute average block-local density de_avg
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
  const enzo_float de_avg = de_sum/count;

  if (de_avg > overdensity_threshold_) {
    // Set mask = 1 everywhere if above threshold
    for (int iz=gz; iz<mz-gz; iz++) {
      int kz=((iz-gz)*nz)/(mz-2*gz);
      for (int iy=gy; iy<my-gy; iy++) {
        int ky=((iy-gy)*ny)/(my-2*gy);
        for (int ix=gx; ix<mx-gx; ix++) {
          int kx=((ix-gx)*nx)/(mx-2*gx);
          int k=kx + nx*(ky + ny*kz);
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

  int level = block->level();

  // Create level array chare elements if this is the base level
  if (level == level_base_) {
    create_level_arrays_(block);
  }

  // Forward count and mask to parent to merge

  if (level == 0) {

    // if root level send final counts to root simulation
    proxy_enzo_simulation[0].p_infer_set_array_count(count);

  } else { // level > 0

    Index index_parent = block->index().index_parent();

    auto parent_block = enzo::block_array()[index_parent];

    if (level > level_base_) {

      // if level > level_base, accumulate level array masks
      int ic3[3] = {0,0,0};
      block->index().child(level,ic3,ic3+1,ic3+2);
      int mx,my,mz;
      std::tie(mx,my,mz) = mask_dims_(block->level());
      const int n = mx*my*mz;
      char ** mask = scalar_mask_(block);

      // merge mask into parent's
      parent_block.p_method_infer_merge_masks (n,*mask,ic3);

    } else { // 0 < level <= level_base_

      // if level <= level base, count inference arrays

      if (level == level_base_) {
        // inference arrays associated with level_base_ blocks
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
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

  // Return control back to EnzoMethodInference
  method->merge_masks(this,n,mask,ic3);
}

//----------------------------------------------------------------------

void EnzoMethodInference::merge_masks
(Block * block, int n_in, char *mask_in, int ic3[3])
{
  const int level = block->level();

  // Merge masks here only if called by a child block
  // (also called by self for synchronization)

  if (n_in > 0) {

    // allocate this block's mask if not allocated yet
    int nx,ny,nz;
    std::tie(nx,ny,nz) = mask_dims_(level);
    mask_allocate_(block,nx,ny,nz);
    char * mask = *scalar_mask_(block);

    if (level < level_array_) {

      // merge multiple masks into one if arrays smaller than block
      merge_masks_2to1_ (mask, mask_in, level, ic3);

    } else { // (level >= level_array)

      // or logical-or masks if blocks smaller than arrays
      merge_masks_1to1_ (mask, mask_in, level);
    }
  }

  // Check if received all child block masks
  if (sync_child_(block).next()) {
    // if so, continue to root level
    forward_create_array_ (block,0);
  }
}

//----------------------------------------------------------------------

void EnzoMethodInference::merge_masks_2to1_
(char * mask, char * mask_in, int level, int ic3[3])
{
  // get my mask size
  int nx,ny,nz;
  std::tie(nx,ny,nz) = mask_dims_(level);

  // get child's mask size
  const int ncx = std::max(nx/2,1);
  const int ncy = std::max(ny/2,1);
  const int ncz = std::max(nz/2,1);

  // copy child's mask into appropriate mask corner
  for (int icz=0; icz<ncz; icz++) {
    int iz = ic3[2]*ncz+icz;
    for (int icy=0; icy<ncy; icy++) {
      int iy = ic3[1]*ncy+icy;
      for (int icx=0; icx<ncx; icx++) {
        int ix = ic3[0]*ncx+icx;

        int i  = ix  + nx *(iy  + ny *iz);
        int ic = icx + ncx*(icy + ncy*icz);

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
  int nx,ny,nz;
  std::tie(nx,ny,nz) = mask_dims_(level);

  // merge child's mask into my mask using logical-or (ie max)
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
        int i = ix + nx*(iy + ny*iz);
        // merge mask values
        mask[i] = std::max(mask[i],mask_in[i]);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_infer_set_array_count(int count)
{
  // Update global count of inference arrays
  infer_count_arrays_ += count;

  // After receiving all expected counts from level-0 blocks...
  if (sync_infer_count_.next()) {

    // initialize synchronization counters,
    // (+1 self-sync in case no inference arrays)
    sync_infer_create_.set_stop(infer_count_arrays_ + 1);
    sync_infer_done_.set_stop(infer_count_arrays_);

    // clear counter for next call,
    infer_count_arrays_ = 0;

    // +1 self-sync call to check_create()
    infer_check_create_();
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_infer_array_created()
// Count number of inference arrays (level array elements) created.
// Called by EnzoLevelArray constructor
{
  infer_check_create_();
}

//----------------------------------------------------------------------
void EnzoSimulation::infer_check_create_()
{
  // Wait for all level array chares to be created...
  if (sync_infer_create_.next()) {

    // If no inference arrays (just +1 call)
    if (sync_infer_create_.stop() == 1) {

      // then exit
      enzo::block_array().p_method_infer_exit();

    } else {

      // otherwise finalize level array chare array
      proxy_level_array.doneInserting();

      // and have level array elements request data from blocks
      proxy_level_array.p_request_data();
    }
    // reset sync counter for next call
    sync_infer_create_.reset();
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::p_request_data()
{
  // Get the index of the (unique) block in level_base_ that overlaps
  // this element
  Index index_block = get_block_index_();

  const int il3[3] = {thisIndex[0],thisIndex[1],thisIndex[2]};

  // Send request for data to block (which may forward the request to
  // its children until leaf blocks reached)
  enzo::block_array()[index_block].p_method_infer_request_data(il3);
}

//----------------------------------------------------------------------

Index EnzoLevelArray::get_block_index_()
{
  Index index;

  // Get "array" part of "array-of-octrees" block index
  int iax = thisIndex[0] >> level_array_;
  int iay = thisIndex[1] >> level_array_;
  int iaz = thisIndex[2] >> level_array_;

  index.set_array(iax,iay,iaz);

  // get "octree" part of "array-of-octrees" block index
  for (int level=1; level<=level_base_; level++) {
    int tree_shift = (level_array_ - level);
    int ic3[3] = {
      (thisIndex[0] >> tree_shift) & 1,
      (thisIndex[1] >> tree_shift) & 1,
      (thisIndex[2] >> tree_shift) & 1};
    index.push_child(ic3[0],ic3[1],ic3[2]);
  }

  // return the Block index corresponding to this EnzoLevelArray
  // element
  return index;
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_request_data (int ia3[3])
{
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

  // Return control back to EnzoMethodInference
  method->request_data(this,ia3);
}

//----------------------------------------------------------------------

void EnzoMethodInference::request_data (Block * block, int ia3[3])
{

  Index index_block = block->index();
  const int level = block->level();

  if (block->is_leaf()) {

    Field field = block->data()->field();

    // determine refinement factors (rx,ry,rz) based on rank
    const int rank = cello::rank();
    const int rx = (rank >= 1) ? 2 : 1;
    const int ry = (rank >= 2) ? 2 : 1;
    const int rz = (rank >= 3) ? 2 : 1;

    // determine "extra" zones required based on rank
    const int ex = (rank >= 1) ? 1 : 0;
    const int ey = (rank >= 2) ? 1 : 0;
    const int ez = (rank >= 3) ? 1 : 0;

    // Verify field_group is not empty
    ASSERT1 ("EnzoMethodInference::request_data()",
             "Field group %s is empty or does not exist",
             field_group_.c_str(),
             (num_fields_ > 0));

    // Serialize data to send to requesting level array element

    // Compute field portion offsets and sizes for /all/ fields
    // (stored in arrays to avoid multiple calls to
    // get_block_portion_() )

    const int n = num_fields_;

    std::vector<int> ox(n),oy(n),oz(n);    // field portion offsets
    std::vector<int> nx(n),ny(n),nz(n);    // field portion sizes
    std::vector<int> nxa(n),nya(n),nza(n);   // array portion sizes
    std::vector<int> oxa(n),oya(n),oza(n);   // array portion sizes

    // compute buffer size nb
    int nb = 0;
    for (int i_f=0; i_f<n; i_f++) {

      const std::string field_name =
        cello::field_groups() -> item(field_group_,i_f);

      const int index_field = field.field_id (field_name);

      std::tie(ox[i_f],oy[i_f],oz[i_f],nx[i_f],ny[i_f],nz[i_f]) =
        get_block_portion_(index_block, index_field, ia3);

      nxa[i_f] = nx[i_f];
      nya[i_f] = ny[i_f];
      nza[i_f] = nz[i_f];
      // determine array offsets if needed
      if (level > level_array_) {
        int cx,cy,cz;
        block->index().child(level,&cx,&cy,&cz);

        unsigned factor = 1 << (level - level_array_);
        unsigned mask = factor - 1;
        oxa[i_f] = (cx & mask)*m3_infer_[0]/factor;
        oya[i_f] = (cy & mask)*m3_infer_[1]/factor;
        oza[i_f] = (cz & mask)*m3_infer_[2]/factor;
      }
      for (int l=level; l>level_infer_; l--) {
        nxa[i_f] = (nxa[i_f]-2*ex)/rx+2*ex;
        nya[i_f] = (nya[i_f]-2*ey)/ry+2*ey;
        nza[i_f] = (nza[i_f]-2*ez)/rz+2*ez;
      }
      // Reserve storage for array size (3) offsets (3), and field
      // value portion (taking into account size after any restrict
      // operations)
      nb += 6 + nxa[i_f]*nya[i_f]*nza[i_f];
    }

    // allocate buffer
    std::vector< enzo_float> buffer_values;
    buffer_values.resize(nb);

    // copy field portions

    if (level > level_infer_) {

      int i_b = 0;

      // for each field in the field group
      for (int i_f=0; i_f<n; i_f++) {

        // Get the field index and valuesfield_name
        const std::string field_name =
          cello::field_groups() -> item(field_group_,i_f);

        enzo_float * field_values = (enzo_float *)field.values(field_name);

        const int index_field = field.field_id (field_name);

        // Find the dimension of the field (may vary between fields
        // depending on centering)
        int mx,my,mz;
        field.dimensions (index_field,&mx,&my,&mz);

        // determine offset into field
        const int ofx = ox[i_f];
        const int ofy = oy[i_f];
        const int ofz = oz[i_f];
        const int of = ofx + mx*(ofy + my*ofz);

        // pointer to field portion to restrict
        enzo_float * a_f = nullptr;
        enzo_float * a_c = nullptr;
        int ncx,ncy,ncz;
        int nfx,nfy,nfz;
        int mfx,mfy,mfz;

        // Apply linear restriction level - level_infer_ times
        for (int l=level; l>level_infer_; l--) {

          // restrict level l to l-1
          const bool is_first = (l == level);
          const bool is_last  = (l == level_infer_ + 1);

          // Input: field if first, else previous output
          if (is_first) {
            a_f = field_values + of;
            nfx = nx[i_f];
            nfy = ny[i_f];
            nfz = nz[i_f];
            mfx = mx;
            mfy = my;
            mfz = mz;
          } else {
            a_f = a_c;
            nfx = (nfx-2*ex)/rx+2*ex;
            nfy = (nfy-2*ey)/ry+2*ey;
            nfz = (nfz-2*ez)/rz+2*ez;
            mfx = nfx;
            mfy = nfy;
            mfz = nfz;
          }

          // Output: array if last, else temporary
          ncx = (nfx-2*ex)/rx+2*ex;
          ncy = (nfy-2*ey)/ry+2*ey;
          ncz = (nfz-2*ez)/rz+2*ez;
          if (is_last) {
            buffer_values[i_b++] = ncx;
            buffer_values[i_b++] = ncy;
            buffer_values[i_b++] = ncz;
            buffer_values[i_b++] = oxa[i_f];
            buffer_values[i_b++] = oya[i_f];
            buffer_values[i_b++] = oza[i_f];
            a_c = buffer_values.data() + i_b;
            i_b += ncx*ncy*ncz;
          } else {
            a_c = new enzo_float [ncx*ncy*ncz];
          }
#ifdef DEBUG_INFER
          CkPrintf ("DEBUG_INFER coarsen a_c %d %d %d  %d %d %d  %d %d %d\n",
                    ncx,ncy,ncz,ncx,ncy,ncz,ex,ey,ez);
          CkPrintf ("DEBUG_INFER coarsen a_f %d %d %d  %d %d %d  %d %d %d\n",
                    mfx,mfy,mfz,nfx,nfy,nfz,ex,ey,ez);
#endif
          coarsen_(a_c,ncx,ncy,ncz,ncx,ncy,ncz,ex,ey,ez,
                   a_f,mfx,mfy,mfz,nfx,nfy,nfz,ex,ey,ez);

#ifdef DEBUG_INFER
        int c,cx,cy,cz;
        block->index().child(block->level(),&cx,&cy,&cz);
        c=1+cx+2*(cy+2*cz);
        for (int iz=0; iz<ncz; iz++) {
          for (int iy=0; iy<ncy; iy++) {
            for (int ix=0; ix<ncx; ix++) {
              const int i = ix + ncx*(iy+ncy*iz);
              a_c[i] = c;
            }
          }
        }
#endif

        // delete temporary when done
          if (!is_first) delete [] a_f;

        }
      }
      ASSERT2("EnzoMethodInference::request_data",
              "Mismatch betwen expected %d and actual %d buffer size",
              i_b,nb,
              (i_b == nb));

    } else {

      // copy fields data

      int i_b = 0; // buffer index

      for (int i_f=0; i_f<n; i_f++) {

        const std::string field_name =
          cello::field_groups() -> item(field_group_,i_f);

        enzo_float * field_values = (enzo_float *)field.values(field_name);

        const int index_field = field.field_id (field_name);

        int mx,my,mz;
        field.dimensions (index_field,&mx,&my,&mz);

        // First copy o[xyz] and n[xyz] indices to buffer
        // (to avoid having to recompute)

        buffer_values[i_b++] = nx[i_f];
        buffer_values[i_b++] = ny[i_f];
        buffer_values[i_b++] = nz[i_f];
        buffer_values[i_b++] = oxa[i_f];
        buffer_values[i_b++] = oya[i_f];
        buffer_values[i_b++] = oza[i_f];

        // Then copy field portion values to buffer
#ifdef DEBUG_INFER
        int c,cx,cy,cz;
        block->index().child(block->level(),&cx,&cy,&cz);
        c=1+cx+2*(cy+2*cz);
#endif
        for (int iz=0; iz<nz[i_f]; iz++) {
          const int ifz = oz[i_f] + iz;
          for (int iy=0; iy<ny[i_f]; iy++) {
            const int ify = oy[i_f] + iy;
            for (int ix=0; ix<nx[i_f]; ix++) {
              const int ifx = ox[i_f] + ix;
              const int iff = ifx + mx*(ify+my*ifz);
#ifdef DEBUG_INFER
              buffer_values[i_b++] = c;
#else
              buffer_values[i_b++] = field_values[iff];
#endif
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

#ifdef TRACE_INFER
    CkPrintf ("TRACE_INFER %s p_transfer_data SEND %d %d %d\n",
              block->name().c_str(),ia3[0],ia3[1],ia3[2]);
#endif
    proxy_level_array[index3].p_transfer_data
      (index_block,nb,buffer_values.data());

  } else { // not leaf

    // forward request to all intersecting child blocks
    ItChild it_child(cello::rank());
    int ic3[3];
    while (it_child.next(ic3)) {
      Index index_child = index_block.index_child(ic3);
      if (block_intersects_array_(index_child,ia3)) {
        enzo::block_array()[index_child].p_method_infer_request_data(ia3);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::p_transfer_data
(Index index, int nb, enzo_float * buffer_values )
{
#ifdef TRACE_INFER
    CkPrintf ("TRACE_INFER p_transfer_data RECV %d %d %d\n",
              thisIndex[0],thisIndex[1],thisIndex[2]);
#endif
  // Loop through fields
  int i_b = 0;
  const int level = index.level();
  const int rank = cello::rank();
  const int rx = (rank >= 1) ? 2 : 1;
  const int ry = (rank >= 2) ? 2 : 1;
  const int rz = (rank >= 3) ? 2 : 1;

  for (int i_f=0; i_f<num_fields_; i_f++) {

    const int nbx = buffer_values[i_b++];
    const int nby = buffer_values[i_b++];
    const int nbz = buffer_values[i_b++];
    const int oax = buffer_values[i_b++];
    const int oay = buffer_values[i_b++];
    const int oaz = buffer_values[i_b++];

    const std::string field_name =
      cello::field_groups() -> item(field_group_,i_f);

    const int index_field = cello::field_descr()->field_id (field_name);

    // Get initial field portion and final array values
    enzo_float * field = buffer_values + i_b;
    enzo_float * array = field_values_[i_f].data();

    enzo_float * a_c = field;
    enzo_float * a_f = array;
    int ecx = (rank >= 1) ? 1 : 0;
    int ecy = (rank >= 2) ? 1 : 0;
    int ecz = (rank >= 3) ? 1 : 0;
    int efx = ecx;
    int efy = ecy;
    int efz = ecz;

    if (level < level_infer_) {
      // Interpolate
      // initialize coarse size
      int ncx,ncy,ncz;
      int nfx,nfy,nfz;
      for (int l=level; l<level_infer_; l++) {

        const bool is_first = (l == level);
        const bool is_last  = (l == level_infer_ - 1);

        // Input: field if first, else previous output
        if (is_first) {
          ncx = nbx;
          ncy = nby;
          ncz = nbz;
          a_c = field;
        } else {
          ncx = (ncx-2*ecx)*rx+2*ecx;
          ncy = (ncy-2*ecy)*ry+2*ecy;
          ncz = (ncz-2*ecz)*rz+2*ecz;
          a_c = a_f;
        }

        // Output: array if last, else temporary
        int mfx,mfy,mfz;
        if (is_last) {
          nfx = (ncx-2*efx)*rx;
          nfy = (ncy-2*efy)*ry;
          nfz = (ncz-2*efz)*rz;
          mfx = nix_;
          mfy = niy_;
          mfz = niz_;
          a_f = array + oax + mfx*(oay + mfy*oaz);
          efx = 0;
          efy = 0;
          efz = 0;
        } else {
          nfx = (ncx-2*efx)*rx+2*efx;
          nfy = (ncy-2*efy)*ry+2*efy;
          nfz = (ncz-2*efz)*rz+2*efz;
          mfx = nfx;
          mfy = nfy;
          mfz = nfz;
          a_f = new enzo_float [mfx*mfy*mfz];
          efx = efx;
          efy = efy;
          efz = efz;
        }

#ifdef DEBUG_INFER
        CkPrintf ("DEBUG_INFER interpolate a_f %d %d %d  %d %d %d  %d %d %d\n",
                  mfx,mfy,mfz,nfx,nfy,nfz,efx,efy,efz);
        CkPrintf ("DEBUG_INFER interpolate a_c %d %d %d  %d %d %d  %d %d %d\n",
                  ncx,ncy,ncz,ncx,ncy,ncz,ecx,ecy,ecz);
#endif
        interpolate_(a_f,mfx,mfy,mfz,nfx,nfy,nfz,efx,efy,efz,
                     a_c,ncx,ncy,ncz,ncx,ncy,ncz,ecx,ecy,ecz);

        // delete temporary when done
        if (!is_first) delete [] a_c;

      }

    } else { // copy

      // Copy (if restriction is required it's performed at source)
      int l = level;
      const int mfx = nix_;
      const int mfy = niy_;
      const int mfz = niz_;
      const int ncx = nbx;
      const int ncy = nby;
      const int ncz = nbz;
      a_f = array + oax + mfx*(oay + mfy*oaz);
#ifdef DEBUG_INFER
      CkPrintf ("DEBUG_INFER copy a_d %d %d %d  %d %d %d  %d %d %d\n",
                mfx,mfy,mfz,nbx-2,nby-2,nbz-2,0,0,0);
      CkPrintf ("DEBUG_INFER copy a_s %d %d %d  %d %d %d  %d %d %d\n",
                ncx,ncy,ncz,ncx-2,ncy-2,ncz-2,0,0,0);
#endif
      copy_(a_f,mfx,mfy,mfz,nbx-2,nby-2,nbz-2,
            a_c,ncx,ncy,ncz,ncx-2,ncy-2,ncz-2);
    }

    // advance pointer to next field
    i_b += nbx*nby*nbz;

  }

  ASSERT2("EnzoLevelArray::p_transfer_data()",
          "Mismatch betwen expected %d and actual %d buffer size",
          i_b,nb,
          (i_b == nb));

  // Check if this is the last incoming block portion expected by
  // tracking volumes

  int level_volume = std::max(0,index.level() - level_array_);
  float volume_portion = pow(1.0/cello::num_children(),level_volume);
  volume_ratio_ += volume_portion;

  ASSERT1 ("EnzoLevelArray::p_transfer_data()",
           "volume_ratio_ %g not between 0.0 and 1.0 as expected",
           volume_ratio_,
           (0.0 <= volume_ratio_) && (volume_ratio_ <= 1.0));

  if (volume_ratio_ >= 1.0) {
    // reset volume counter for next call
    volume_ratio_ = 0.0;
    // When done, call inference (note safe to compare float with constant
    // since guaranteed no roundoff error)
    apply_inference();
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::apply_inference()
{
  double lower[3],upper[3];
  this->lower(lower);
  this->upper(upper);

#ifdef DEBUG_INFER
  CkPrintf ("DEBUG_INFER level array count na[xyz] %d %d %d\n",nax_,nay_,naz_);
  CkPrintf ("DEBUG_INFER level array size  ni[xyz] %d %d %d\n",nix_,niy_,niz_);
  for (int i_f=0; i_f<num_fields_; i_f++) {
    std::string field_name=cello::field_groups() -> item(field_group_,i_f);
    double sum=0.0;
    double min=+1e10;
    double max=-1e10;
    double wsum=0.0;

    for (int iz=0; iz<niz_; iz++) {
      for (int iy=0; iy<niy_; iy++) {
        for (int ix=0; ix<nix_; ix++) {
          const int i = ix+nix_*(iy+niy_*iz);
          sum += field_values_[i_f][i];
          min = std::min(min,field_values_[i_f][i]);
          max = std::max(max,field_values_[i_f][i]);
          wsum += (ix+1+2*(iy+1)-3*(iz+1))*field_values_[i_f][i];
        }
      }
    }
    CkPrintf ("DEBUG_INFER %s SUM %d %d %d %g\n",
              field_name.c_str(),thisIndex[0],thisIndex[1],thisIndex[2],sum);
    CkPrintf ("DEBUG_INFER %s MIN %d %d %d %g\n",
              field_name.c_str(),thisIndex[0],thisIndex[1],thisIndex[2],min);
    CkPrintf ("DEBUG_INFER %s MAX %d %d %d %g\n",
              field_name.c_str(),thisIndex[0],thisIndex[1],thisIndex[2],max);
    CkPrintf ("DEBUG_INFER %s WSUM %d %d %d %g\n",
              field_name.c_str(),thisIndex[0],thisIndex[1],thisIndex[2],wsum);
    fflush(stdout);
  }
#endif
  //==================================================
  //
  // ADD DEEP LEARNING INFERENCE HERE
  //
  //==================================================

  const EnzoConfig * enzo_config = enzo::config();

// if inference_method == "starnet"
// TUTORIAL: https://github.com/prabhuomkar/pytorch-cpp/blob/master/tutorials/basics/pytorch_basics/main.cpp
    
  double sphere_x=0.0, sphere_y=0.0, sphere_z=0.0;
  double sphere_r=0.0;
  double yield_SNe=0.0, yield_HNe=0.0, yield_PISNe=0.0;

  #ifdef CONFIG_USE_TORCH
    // StarFind:
    // ----------------------
    // prepare fields, shove into torch::Tensor "sample" with dimensions {1, N_fields, dim, dim, dim} 
    // loop through fields in "Inference" group
    //   see StarNetDataLoader()
    // density, H2I_density, total_energy, div(v), and metal_density

    const std::string stage1_checkpoint = "/home1/07320/whick002/StarNetRuntime/model_checkpoints/smalldense.jtpt";
    
    // Deserialized stage-1 checkpoint (U-net)
    torch::jit::script::Module stage1 = starfind_->stage1(); 
    torch::Tensor field_data = torch::zeros({1,5,nix_,niy_,niz_}, torch::dtype(torch::kFloat32)); // array that gets passed into model

    // array increments
    const int idx = 1;
    const int idy = nix_;
    const int idz = nix_*niy_;

    // create velocity divergence derived field
    enzo_float * velocity_divergence = new enzo_float[nix_*niy_*niz_];
    enzo_float * vx = field_values_[2].data();
    enzo_float * vy = field_values_[3].data();
    enzo_float * vz = field_values_[4].data();

    int nx,ny,nz;
    cello::hierarchy()->root_size(&nx,&ny,&nz);

    const double hx = (upper[0] - lower[0]) / nix_;
    const double hy = (upper[1] - lower[1]) / niy_;
    const double hz = (upper[2] - lower[2]) / niz_;

    for (int iz=0; iz<niz_; iz++){
      for (int iy=0; iy<niy_; iy++){
        for (int ix=0; ix<nix_; ix++){
          int i = INDEX(ix,iy,iz,nix_,niy_);

          double divVx=0.0, divVy=0.0, divVz=0.0;

          if (ix == 0)         divVx = vx[i+idx]-vx[i];
          else if (ix == nix_-1) divVx = vx[i]-vx[i-idx];
          else                 divVx = 0.5 * (vx[i+idx] - vx[i-idx]);

          if (iy == 0)         divVy = vy[i+idy]-vy[i];
          else if (iy == niy_-1) divVy = vy[i]-vy[i-idy];
          else                 divVy = 0.5 * (vy[i+idy] - vy[i-idy]);

          if (iz == 0)         divVz = vz[i+idz]-vz[i];
          else if (iz == niz_-1) divVz = vz[i]-vz[i-idz];
          else                 divVz = 0.5 * (vz[i+idz] - vz[i-idz]);

          velocity_divergence[i] = divVx + divVy + divVz;

        } // ix
      } // iy
    } // iz
 
    // scale field such that x = (x-mean) / std
    EnzoUnits * enzo_units = enzo::units();
    double rhounit = enzo_units->density();
    double vunit = enzo_units->velocity();
    double lunit = enzo_units->length();
    double Zunit = enzo_constants::metallicity_solar;
 
    // means and stds from StarNetRuntime/resources/Scaling_64x160.torch (generated from the Phoenix simulations) 
    // All values are in CGS units except for total_energy
    std::vector <double> means = {3.395644060629469e-27/rhounit, 1.5020361582642288e-30/rhounit, 
                -793.0158410664843/vunit, 
                96878457078.45897, 2.3778995827734495e-08, 2.419114340099057e-07};
    std::vector <double> stds  = {4.858201188227477e-26/rhounit, 7.886306542228033e-29/rhounit, 
                8358.282154906454/vunit,
                531795315558.7779, 2.5806363586099314e-06, 2.9415985082639072e-05};

#ifdef PRINT_FIELD_VALUES
  for (int i_f = 0; i_f < num_fields_; i_f++) {
      // print out all field values
      for (int iz=0; iz<niz_; iz++){
        for (int iy=0; iy<niy_; iy++){
          for (int ix=0; ix<nix_; ix++){
            int i = ix + nix_*(iy + niy_*iz);
            CkPrintf("field_%d(%d,%d,%d) = %e\n",i_f,ix,iy,iz, field_values_[i_f].data()[i]);
          }
        }
      }
  }
#endif

    // i_f = iterator for fields in "inference" group
    // i_f_ = iterator for fields going into model
    int i_f_ = 0; 
    for (int i_f=0; i_f < num_fields_-2; i_f++) { 
      // num_fields-2 because PopIII_metal_density is only used for scaling the metal field -- it 
      // doesn't actually get passed into the model directly. Will handle this case separately

      enzo_float * array;
      if (i_f == 2) {
        array = velocity_divergence; // vx, vy, and vz are all included in this, so skip i=2-4
        i_f = 4; // i_f = 6 is metal_density
      } 

      else {
          array = field_values_[i_f].data();
      }

      for (int iz=0; iz<niz_; iz++){
        for (int iy=0; iy<niy_; iy++){
          for (int ix=0; ix<nix_; ix++){
            int i = ix + nix_*(iy + niy_*iz);
            // see https://pytorch.org/cppdocs/notes/tensor_indexing.html
            field_data.index_put_( {0,i_f_,ix,iy,iz}, 
                         ((double) array[i] - means[i_f_])/stds[i_f_] );

          } // ix
        } // iy
      } // iz
      i_f_ += 1;
    } // i_field

    // Now deal with metal fields
    enzo_float * array_tot    = field_values_[6].data(); // metal_density
    enzo_float * array_popIII = field_values_[7].data(); // PopIII_metal_density
    double mean_popII  = means[4], std_popII  = stds[4];
    double mean_popIII = means[5], std_popIII = stds[5];
    for (int iz=0; iz<niz_; iz++){
      for (int iy=0; iy<niy_; iy++){
        for (int ix=0; ix<nix_; ix++){
          int i = ix + nix_*(iy + niy_*iz);
          double rho_popIII = std::min(array_popIII[i], array_tot[i]); // rho_popIII <= rho_total always
          double rho_popII  = array_tot[i] - array_popIII[i];

          double val = 0.0;
          // see StarNetRuntime/StarNetDataLoader.py
          val  = (rho_popIII-mean_popIII)*std_popII;
          //val -= mean_popIII*std_popII;
          val += (rho_popII-mean_popII)*std_popIII;
          val /= std_popIII*std_popII;

          //CkPrintf("scaled metal density: %1.2e; metal_density: %1.2e; rho_popIII: %1.2e\n", val, array_tot[i], rho_popIII);

          // see https://pytorch.org/cppdocs/notes/tensor_indexing.html
          field_data.index_put_( {0,4,ix,iy,iz}, val );
        }
      }
    }

    delete [] velocity_divergence;

    std::vector<torch::jit::IValue> sample; 
    sample.push_back(field_data);

    // if skipping stage 1 classifier, all inference blocks
    // are treated as potential star-forming blocks 
    // and will be fed into the stage 2 inception network
    bool star_in_block = true;

    if (enzo_config->method_inference_starnet_S1) {
      // see step 4 of https://pytorch.org/tutorials/advanced/cpp_export.html 

      at::Tensor prediction_s1 = stage1.forward(sample).toTensor(); //(pstar, pnostar)
    
      double num_0 = std::exp(prediction_s1.index({0,0}).item<double>());
      double num_1 = std::exp(prediction_s1.index({0,1}).item<double>());
      double pnostar = num_0 / (num_0 + num_1);
      double pstar   = num_1 / (num_0 + num_1);
 
      CkPrintf("EnzoMethodInference::S1 prediction = (%f, %f), (pnostar, pstar) = (%f, %f); pnostar + pstar = %f\n", prediction_s1.index({0,0}).item<double>(), prediction_s1.index({0,1}).item<double>(), pnostar, pstar, pstar+pnostar);

      //ASSERT("EnzoMethodInference::apply_inference()",
      //       "pstar+pnostar != 1.0", pstar+pnostar == 1.0);
        
      star_in_block = (pstar > pnostar);
    }

    bool SF_confirmed = false; // flag for if S2 network identified star-forming regions

    if (star_in_block && enzo_config->method_inference_starnet_S2) {

      CkPrintf("EnzoMethodInference::Calling S2 inception U-net...\n");
      
      // Deserialized stage-2 checkpoint (inception net)
      torch::jit::script::Module stage2 = starfind_->stage2(); // shape = {1,2,nx,ny,nz}

      at::Tensor prediction_s2 = stage2.forward(sample).toTensor();

      int N_edge = 0; // number of edge layers to consider when evaluating S2 predictions
      double pstar_i, pnostar_i, num_0_i, num_1_i;
      double rho_x = 0.0, rho_y = 0.0, rho_z = 0.0, rho = 0.0; // for calculating mean position

      enzo_float * density = field_values_[0].data();
      enzo_float * metal_density = field_values_[6].data();

      bool is_valid = true;
      double zcrit = 3.1e-6; // TODO: Make this a parameter!!!!!

      for (int iz=0; iz<niz_; iz++){
       for (int iy=0; iy<niy_; iy++){
         for (int ix=0; ix<nix_; ix++){
           int i = INDEX(ix,iy,iz,nix_,niy_);

           num_0_i = std::exp(prediction_s2.index({0,0,ix,iy,iz}).item<double>());
           num_1_i = std::exp(prediction_s2.index({0,1,ix,iy,iz}).item<double>());
           pnostar_i = num_0_i / (num_0_i + num_1_i);
           pstar_i   = num_1_i / (num_0_i + num_1_i);

           if (pstar_i > pnostar_i) {
             SF_confirmed = true;

             double x_i = lower[0]+(ix+0.5)*hx;
             double y_i = lower[1]+(iy+0.5)*hy;
             double z_i = lower[2]+(iz+0.5)*hz;

             // contribute to the mean predicted position
             double rho_i = density[i];
             rho += rho_i;
             rho_x += rho_i * x_i;
             rho_y += rho_i * y_i;
             rho_z += rho_i * z_i;

             if (metal_density[i]/rho_i / enzo_constants::metallicity_solar > zcrit) {
               // if a flagged cell is above the critical metallicity separating PopII from PopIII SF,
               // this identified region is invalid. Don't do anything in this case
               is_valid = false;
             } 

             // if a flagged cell is within N_edge of the border of this inference block,
             // this identified region is also invalid.
             
             is_valid *= (N_edge <= iz) && (iz < niz_-N_edge);
             is_valid *= (N_edge <= iy) && (iy < niy_-N_edge);
             is_valid *= (N_edge <= ix) && (ix < nix_-N_edge);

             CkPrintf("EnzoMethodInference::Starfind predicts potential Pop III star formation at (x, y, z) = (%f,%f,%f); index = %d (%d,%d,%d); lower = (%f,%f,%f); upper = (%f,%f,%f)\n", x_i, y_i, z_i, i, ix, iy, iz, lower[0], lower[1], lower[2], upper[0], upper[1], upper[2]);
           } 
         } // ix
       } // iy
      } // iz

      // REMOVE THIS LINE!!!!!! JUST FOR DEBUGGING
      //is_valid = true;

      bool method_inference_starnet_feedback = true; // TODO: Make this a parameter
      // only do FB if S2 network identifies star-forming regions
      if ( (SF_confirmed*is_valid) && method_inference_starnet_feedback) {
        // center of FB sphere is (xmean, ymean, zmean)
        double rho_inv = 1.0 / rho;
        sphere_x = rho_x * rho_inv;
        sphere_y = rho_y * rho_inv;
        sphere_z = rho_z * rho_inv;
     
        // determine properties of population and metal yields
        std::vector<double> masses, creationtimes;
      
        int Nstars = fbnet_->get_Nstars();
        for (int i=0; i<Nstars; i++) {
          double mass = fbnet_->get_mass();
          if ((11.0 <= mass) && (mass < 20.0)) {
            yield_SNe += fbnet_->metal_yield_SNe(mass);
          }
          else if ((20.0 <= mass) && (mass < 40.0)) {
            yield_HNe += fbnet_->metal_yield_HNe(mass);
          }
          else if ((140.0 <= mass) && (mass < 260.0)) {
            yield_PISNe += fbnet_->metal_yield_PISNe(mass);
          }
          masses.push_back(mass);
          creationtimes.push_back(fbnet_->get_creationtime());
        }
 
        // The evolved remnant radii in the Phoenix simulations were determined
        // by expanding a sphere until the average metallicity within
        // drops below Z_crit. Since the evolved remnants won't be perfectly
        // spherical in general, and since metallicities could in principle
        // be much larger than Z_crit, the radii returned by the model
        // will be too big. Correct for this by introducing a "radius_modifier" factor.
        // This is introduced in FBNet.py (apply_spherical_feedback). Azton's final
        // simulation that is comparable to the Phoenix simulations used a value of 0.2
        double radius_modifier = 0.2; // TODO: Make this a parameter!
        // calculate radius of evolved remnant
        sphere_r = (yield_SNe+yield_PISNe+yield_HNe > 0) ? 
                 radius_modifier*fbnet_->get_radius( masses, creationtimes ) : 0.0;
        CkPrintf("EnzoMethodInference::FBNet predicts a region with radius %1.2e kpc and metal yields (%1.2e, %1.2e, %1.2e) Msun\n", sphere_r, yield_SNe, yield_HNe, yield_PISNe);
        // put everything into code units
        EnzoUnits * enzo_units = enzo::units();
        double lunit = enzo_units->length();
        double munit = enzo_units->mass();

        sphere_r *= enzo_constants::kpc_cm / lunit;
        yield_SNe *= enzo_constants::mass_solar / munit;
        yield_HNe *= enzo_constants::mass_solar / munit;
        yield_PISNe *= enzo_constants::mass_solar / munit;
    }
   
  } // if pstar > pnostar in a block

  #endif

  // Update blocks with inference results

  //    center of sphere
  double center[3] = {sphere_x, sphere_y, sphere_z};

  //    put a sphere object there and add it to a list sphere_list
  EnzoObjectFeedbackSphere sphere(center, sphere_r, yield_SNe, yield_HNe, yield_PISNe);
 
  //    allocate buffer
  int n = 0;
  SIZE_SCALAR_TYPE(n,EnzoObjectFeedbackSphere,sphere);
  //    initialize buffer
  char * buffer = new char [n];
  char * pc = buffer;
  SAVE_SCALAR_TYPE(pc,EnzoObjectFeedbackSphere,sphere);

  //    Send data to leaf blocks via base-level block
  Index index_block = get_block_index_();
  const int il3[3] = {thisIndex[0],thisIndex[1],thisIndex[2]};
  enzo::block_array()[index_block].p_method_infer_update(n,buffer,il3);

#ifdef TRACE_INFER
  CkPrintf ("TRACE_INFER rectangle %d %g %g %g %g %g %g\n",
            cello::simulation()->cycle(),
            lower[0],lower[1],lower[2],
            upper[0],upper[1],upper[2]);
#endif
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_infer_update(int n, char * buffer, int il3[3])
{
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

  method->update(this, n, buffer, il3);
}

//----------------------------------------------------------------------

void EnzoMethodInference::update ( Block * block, int n, char * buffer, int il3[3])
{
  // Unpack buffer into sphere object
  EnzoObjectFeedbackSphere sphere;
  char *pc = buffer;
  LOAD_SCALAR_TYPE(pc,EnzoObjectFeedbackSphere,sphere);

  Index index_block = block->index();
  const int level = block->level();

  if (block->is_leaf()) {

    // if inference_method == "starnet"
    // Note that there is a check within spawn_remnant_particle
    // that particle is located inside this leaf block. Necessary
    // because EnzoMethodInference::update() forwards data from
    // the root block to ALL leaves contained inside its volume,
    // so the sphere particle we would like to spawn may not
    // necessarily lie in this specific leaf.
    spawn_remnant_particle(block, &sphere);

    // if leaf block, we're done, tell level array element
    Index3 index3(il3[0],il3[1],il3[2]);

    proxy_level_array[index3].p_done (index_block);

  } else {

    // else if non-leaf block, forward data to intersecting child blocks
    ItChild it_child(cello::rank());
    int ic3[3];
    while (it_child.next(ic3)) {
      Index index_child = index_block.index_child(ic3);
      if (block_intersects_array_(index_child,il3)) {
        enzo::block_array()[index_child].p_method_infer_update(n,buffer,il3);
      }
    }
  }

}

//-----------------------------------------

void EnzoMethodInference::spawn_remnant_particle ( Block * block, 
                                    EnzoObjectFeedbackSphere * sphere)
{
  if (sphere->r() == 0.0) {
    // return if sphere object is empty
    return;
  } 
 
  EnzoBlock * enzo_block = enzo::block(block);

  Field field = enzo_block->data()->field();

  int mx,my,mz;
  field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  int gx,gy,gz;
  field.ghost_depth(0,&gx, &gy, &gz);

  double xm,ym,zm;
  double xp,yp,zp;
  double hx,hy,hz;
  enzo_block->data()->lower(&xm,&ym,&zm);
  enzo_block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);
   
  // get 3D grid index for sphere center - account for ghost zones!!
  int ix = (int) std::floor((sphere->pos(0) - xm) / hx) + gx;
  int iy = (int) std::floor((sphere->pos(1) - ym) / hy) + gy;
  int iz = (int) std::floor((sphere->pos(2) - zm) / hz) + gz;

  // make sure that particle is actually within block
  bool in_block = true;
  in_block = in_block && ((gx <= ix) && (ix <= mx-gx));
  in_block = in_block && ((gy <= iy) && (iy <= my-gy));
  in_block = in_block && ((gz <= iz) && (iz <= mz-gz));
  
  if (! in_block ) {
    return;
  }

  EnzoUnits * enzo_units = enzo::units();
  double lunit = enzo_units->length();
  double munit = enzo_units->mass();

  Particle particle = enzo_block->data()->particle();

  // declare particle arrays
  const int it = particle.type_index("popIIIremnant");

  const int ia_x = particle.attribute_index (it, "x");
  const int ia_y = particle.attribute_index (it, "y");
  const int ia_z = particle.attribute_index (it, "z");
  const int ia_r = particle.attribute_index (it, "r");
  const int ia_vx = particle.attribute_index (it, "vx");
  const int ia_vy = particle.attribute_index (it, "vy");
  const int ia_vz = particle.attribute_index (it, "vz");
  
  const int ia_mSNe = particle.attribute_index(it, "yield_SNe");
  const int ia_mHNe = particle.attribute_index(it, "yield_HNe");
  const int ia_mPISNe = particle.attribute_index(it, "yield_PISNe");

  const int ia_t = particle.attribute_index(it, "creation_time");

  int ib  = 0; // batch counter
  int ipp = 0; // counter

  // insert a particle
  int my_particle = particle.insert_particles(it, 1);

  // get batch number and particle index
  particle.index(my_particle, &ib, &ipp);

  // pointers for particle attribute arrays
  enzo_float * px   = (enzo_float *) particle.attribute_array(it, ia_x, ib);
  enzo_float * py   = (enzo_float *) particle.attribute_array(it, ia_y, ib);
  enzo_float * pz   = (enzo_float *) particle.attribute_array(it, ia_z, ib);
  enzo_float * pr   = (enzo_float *) particle.attribute_array(it, ia_r, ib);
  enzo_float * pvx  = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
  enzo_float * pvy  = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
  enzo_float * pvz  = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
  enzo_float * pmSNe = (enzo_float *) particle.attribute_array(it, ia_mSNe, ib);
  enzo_float * pmHNe = (enzo_float *) particle.attribute_array(it, ia_mHNe, ib);
  enzo_float * pmPISNe = (enzo_float *) particle.attribute_array(it, ia_mPISNe, ib);
  enzo_float * pform  = (enzo_float *) particle.attribute_array(it, ia_t, ib);

  int io = ipp;

  px[io] = sphere->pos(0); // code length
  py[io] = sphere->pos(1);
  pz[io] = sphere->pos(2);
  pr[io] = sphere->r();

  pmSNe[io] = sphere->metal_mass_SNe(); // code mass
  pmHNe[io] = sphere->metal_mass_HNe();
  pmPISNe[io] = sphere->metal_mass_PISNe();

  enzo_float * density = (enzo_float *) field.values("density");
  enzo_float * velocity_x = (enzo_float *) field.values("velocity_x");
  enzo_float * velocity_y = (enzo_float *) field.values("velocity_y");
  enzo_float * velocity_z = (enzo_float *) field.values("velocity_z");

  // average particle velocity over many cells to prevent runaway
  double rhosum = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
  for (int iz_=std::max(0,iz-3); iz_ <= std::min(iz+3,mx); iz_++) {
    for (int iy_=std::max(0,iy-3); iy_ <= std::min(iy+3,my); iy_++) {
      for (int ix_=std::max(0,ix-3); ix_ <= std::min(ix+3,mz); ix_++) {
        int i_ = INDEX(ix_,iy_,iz_,mx,my);
        vx += velocity_x[i_]*density[i_];
        vy += velocity_y[i_]*density[i_];
        vz += velocity_z[i_]*density[i_];
        rhosum += density[i_];
      }
    } 
  }  
  vx /= rhosum;
  vy /= rhosum;
  vz /= rhosum;

  pvx[io] = vx;
  pvy[io] = vy;
  pvz[io] = vz;

  CkPrintf("Spawning remnant particle!!! creation_time = %f \n", block->time());
  pform[io] = block->time();
}

//----------------------------------------------------------------------

void EnzoLevelArray::p_done(Index index)
{
  // add volume of block to volume counter
  int level_volume = std::max(0,index.level() - level_array_);
  float volume_portion = pow(1.0/cello::num_children(),level_volume);
  volume_ratio_ += volume_portion;

  ASSERT1 ("EnzoLevelArray::p_done()",
           "volume_ratio_ %g not between 0.0 and 1.0 as expected",
           volume_ratio_,
           (0.0 <= volume_ratio_) && (volume_ratio_ <= 1.0));

  if (volume_ratio_ >= 1.0) {
    // reset volume counter for next call
    volume_ratio_ = 0.0;
    // tell root Simulation object this level array is done
    proxy_enzo_simulation[0].p_infer_done();
  }
}
//----------------------------------------------------------------------

void EnzoSimulation::p_infer_done()
{
  // count level arrays that are done
  if (sync_infer_done_.next()) {
    // when all are done, reset sync counters,
    sync_infer_done_.reset();
    // ... delete chare array,
    proxy_level_array.ckDestroy();
    // ... and exit method
    enzo::block_array().p_method_infer_exit();
  }
}
//----------------------------------------------------------------------

std::tuple<int,int,int,int,int,int>
EnzoMethodInference::get_block_portion_
(Index index, int index_field, int ia3[3])
{
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

  const int ratio = ratio_infer;

  nbx = (rank >= 1) ? nx / ratio + 2*ex : 1;
  nby = (rank >= 2) ? ny / ratio + 2*ey : 1;
  nbz = (rank >= 3) ? nz / ratio + 2*ez : 1;

  // compute field portion offset
  int obx,oby,obz;

  int gx,gy,gz;
  cello::field_descr()->ghost_depth(index_field,&gx,&gy,&gz);

  obx = kbx*nx/ratio + gx - ex;
  oby = kby*ny/ratio + gy - ey;
  obz = kbz*nz/ratio + gz - ez;

  return {obx,oby,obz,nbx,nby,nbz};
}

//----------------------------------------------------------------------

void EnzoMethodInference::coarsen_
(enzo_float * a_c,
 int mcx, int mcy, int mcz,
 int ncx, int ncy, int ncz,
 int ecx, int ecy, int ecz,
 const enzo_float * a_f,
 int mfx, int mfy, int mfz,
 int nfx, int nfy, int nfz,
 int efx, int efy, int efz)
{
  // initialize index offsets, d[0]=0, d[1] = dx, d[2]=dy, d[3]=dx+dy,
  // etc.
  const int dx = 1;
  const int dy = mfx;
  const int dz = mfx*mfy;
  int d[8];
  for (int iz=0; iz<2; iz++) {
    for (int iy=0; iy<2; iy++) {
      for (int ix=0; ix<2; ix++) {
        const int i=ix+2*(iy+2*iz);
        d[i] = ix*dx + iy*dy + iz*dz;
      }
    }
  }

  // Compute restriction based on problem rank
  const int rank = cello::rank();
  if (rank == 1) {
    for (int icx=0; icx<mcx; icx++) {
      int i_f=icx*2;
      int i_c=icx;
      a_c[i_c] = 0.5*(a_f[i_f+d[0]] +
                      a_f[i_f+d[1]]);
    }
  } else if (rank == 2) {
    for (int icy=0; icy<mcy; icy++) {
      for (int icx=0; icx<mcx; icx++) {
        int i_f=icx*2 + mfx*(icy*2);
        int i_c=icx   + mcx*(icy  );
        a_c[i_c] = 0.25*(a_f[i_f+d[0]] +
                         a_f[i_f+d[1]] +
                         a_f[i_f+d[2]] +
                         a_f[i_f+d[3]]);
      }
    }
  } else if (rank == 3) {
    for (int icz=0; icz<mcz; icz++) {
      for (int icy=0; icy<mcy; icy++) {
        for (int icx=0; icx<mcx; icx++) {
          int i_f=icx*2 + mfx*(icy*2 + mfy*icz*2);
          int i_c=icx   + mcx*(icy   + mcy*icz);
          a_c[i_c] = 0.125*(a_f[i_f+d[0]] +
                            a_f[i_f+d[1]] +
                            a_f[i_f+d[2]] +
                            a_f[i_f+d[3]] +
                            a_f[i_f+d[4]] +
                            a_f[i_f+d[5]] +
                            a_f[i_f+d[6]] +
                            a_f[i_f+d[7]]);
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::interpolate_
(enzo_float * a_f,
 int mfx, int mfy, int mfz,
 int nfx, int nfy, int nfz,
 int efx, int efy, int efz,
 const enzo_float * a_c,
 int mcx, int mcy, int mcz,
 int ncx, int ncy, int ncz,
 int ecx, int ecy, int ecz)
{
  // initialize index offsets, d[0]=0, d[1] = dx, d[2]=dy, d[3]=dx+dy,
  // etc.
  const int dx = 1;
  const int dy = mcx;
  const int dz = mcx*mcy;
  int d[8];
  for (int iz=0; iz<2; iz++) {
    for (int iy=0; iy<2; iy++) {
      for (int ix=0; ix<2; ix++) {
        const int i=ix+2*(iy+2*iz);
        d[i] = ix*dx + iy*dy + iz*dz;
      }
    }
  }

  const enzo_float w[2] = { 0.75, 0.25 };

  if (cello::rank() == 1) {
    for (auto ifx=0; ifx<nfx; ifx++) {
      auto iex = ifx+(1-efx);
      auto icx = iex/2;
      auto wx0 = w[iex&1];
      auto wx1 = 1.0-wx0;
      auto i_f = ifx;
      auto i_c = icx;
      a_f[i_f] = (wx0*a_c[i_c+d[0]] +
                  wx1*a_c[i_c+d[1]]);
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
        a_f[i_f] = (wy0*wx0*a_c[i_c+d[0]] +
                    wy0*wx1*a_c[i_c+d[1]] +
                    wy1*wx0*a_c[i_c+d[2]] +
                    wy1*wx1*a_c[i_c+d[3]]);
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
          a_f[i_f] = (wz0*wy0*wx0*a_c[i_c+d[0]] +
                      wz0*wy0*wx1*a_c[i_c+d[1]] +
                      wz0*wy1*wx0*a_c[i_c+d[2]] +
                      wz0*wy1*wx1*a_c[i_c+d[3]] +
                      wz1*wy0*wx0*a_c[i_c+d[4]] +
                      wz1*wy0*wx1*a_c[i_c+d[5]] +
                      wz1*wy1*wx0*a_c[i_c+d[6]] +
                      wz1*wy1*wx1*a_c[i_c+d[7]]);
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoLevelArray::copy_
(enzo_float * a_dst,
 int mdx, int mdy, int mdz, int ndx,int ndy, int ndz,
 const enzo_float * a_src,
 int msx, int msy, int msz, int nsx,int nsy, int nsz)
{
  ASSERT6("EnzoLevelArray::copy_",
          "Mismatch between expected sizes: src %d %d %d dst %d %d %d\n",
          nsx,nsy,nsz,ndx,ndy,ndz,
          (nsx==ndx) && (nsy==ndy) && (nsz==ndz));
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
  EnzoMethodInference * method =
    static_cast<EnzoMethodInference*> (this->method());

  // Return control back to EnzoMethodInference
  method->count_arrays(this,count);
}

void EnzoMethodInference::count_arrays (Block * block, int count)
{
  scalar_count_(block) += count;

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
  compute_done();
}

//======================================================================

void EnzoMethodInference::create_level_arrays_ (Block * block)
{
  const int level = block->level();

  ASSERT2 ("EnzoMethodInference::create_level_arrays_",
           "Block %s not in expected refinement level %d",
           block->name().c_str(),level,
           level == level_base_);

  int mx,my,mz;
  std::tie(mx,my,mz) = mask_dims_(level);
  int m = mx*my*mz;
  char * mask = *scalar_mask_(block);
  int nax,nay,naz;
  level_array_dims_(&nax,&nay,&naz);
  if (mask != nullptr) {
    int i3[3];
    block->index().index_level (i3, level);
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
  const float r = pow(2,(level_array_ - index.level()));
  int lx = (ib3[0] <= ia3[0] && ia3[0] < ib3[0]+r);
  int ly = (ib3[1] <= ia3[1] && ia3[1] < ib3[1]+r);
  int lz = (ib3[2] <= ia3[2] && ia3[2] < ib3[2]+r);
  bool intersects = lx && ly && lz;
  return intersects;
}

//================================================================
