// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-10
/// @brief    Implementation of the Patch class

#include "cello.hpp"

#include "mesh.hpp"
#include "main.hpp"

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
extern CProxy_SimulationCharm  proxy_simulation;
#endif

//----------------------------------------------------------------------

Patch::Patch
(
#ifndef CONFIG_USE_CHARM
 const Factory * factory,
 const FieldDescr * field_descr,
#endif
 int nx,  int ny,  int nz,
 int nx0, int ny0, int nz0,
 int nbx, int nby, int nbz,
 double xm, double ym, double zm,
 double xp, double yp, double zp,
 int id,
 bool allocate_blocks,
 int process_first, int process_last_plus
) throw()
  :
  id_(id),
#ifdef CONFIG_USE_CHARM
  block_array_(0),
  block_exists_(false),
  block_counter_(0),
#else
  factory_(factory),
#endif
  group_process_(GroupProcess::create(process_first,process_last_plus)),
  layout_ (nbx,nby,nbz)
{

 // Check 
  DEBUG2("zp = %f ID = %d",zp,id_);
  DEBUG2("Patch() %p id = %d",this,this->id());
  if ( ! ((nx >= nbx) && (ny >= nby) && (nz >= nbz))) {
	     
    ERROR6("Patch::Patch", 
	   "Patch size (%d,%d,%d) must be larger than blocking (%d,%d,%d)",
	   nx,ny,nz,nbx,nby,nbz);
  }

  // set layout process range

  layout_.set_process_range(0,group_process_->size());

  size_[0] = nx;
  size_[1] = ny;
  size_[2] = nz;
  DEBUG1 ("patch = %p",this);
  DEBUG1 ("patch size = %d",size_[0]);

  offset_[0] = nx0;
  offset_[1] = ny0;
  offset_[2] = nz0;

  blocking_[0] = nbx;
  blocking_[1] = nby;
  blocking_[2] = nbz;

  lower_[0] = xm;
  lower_[1] = ym;
  lower_[2] = zm;

  upper_[0] = xp;
  upper_[1] = yp;
  upper_[2] = zp;

#ifdef CONFIG_USE_CHARM
  DEBUG2("patch %p ID = %d",this,id_);
  allocate_array_(allocate_blocks);
#else
  allocate_array_(allocate_blocks,field_descr);
#endif

#ifdef CONFIG_USE_CHARM
  proxy_simulation.s_initialize();
#endif
}

//----------------------------------------------------------------------

Patch::~Patch() throw()
{
  deallocate_blocks();
  delete group_process_; group_process_ = 0;
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Patch::s_block(CkCallback callback)
{
  TRACE1("Patch::s_block(%d)",block_counter_.count_curr());

  if (block_counter_.remaining() == 0) {

    callback.send();

  }
}

#endif

//----------------------------------------------------------------------

void Patch::size (int * npx, int * npy, int * npz) const throw()
{
  if (npx) (*npx) = size_[0];
  if (npy) (*npy) = size_[1];
  if (npz) (*npz) = size_[2];
}

//----------------------------------------------------------------------

void Patch::offset (int * nx0, int * ny0, int * nz0) const throw()
{
  if (nx0) (*nx0) = offset_[0];
  if (ny0) (*ny0) = offset_[1];
  if (nz0) (*nz0) = offset_[2];
}

//----------------------------------------------------------------------

void Patch::blocking (int * nbx, int * nby, int * nbz) const throw()
{
  if (nbx) (*nbx) = blocking_[0];
  if (nby) (*nby) = blocking_[1];
  if (nbz) (*nbz) = blocking_[2];
}

//----------------------------------------------------------------------

Layout * Patch::layout () throw()
{
  return &layout_;
}

//----------------------------------------------------------------------

const Layout * Patch::layout () const throw()
{
  return &layout_;
}

//----------------------------------------------------------------------
  
void Patch::lower(double * xm, double * ym, double * zm) const throw ()
{
  if (xm) (*xm) = lower_[0];
  if (ym) (*ym) = lower_[1];
  if (zm) (*zm) = lower_[2];
}

//----------------------------------------------------------------------
void Patch::upper(double * xp, double * yp, double * zp) const throw ()
{
  if (xp) (*xp) = upper_[0];
  if (yp) (*yp) = upper_[1];
  if (zp) (*zp) = upper_[2];
}

//----------------------------------------------------------------------

void Patch::deallocate_blocks() throw()
{

#ifdef CONFIG_USE_CHARM

  if (block_exists_) {
    block_array_->ckDestroy();
    delete block_array_; block_array_ = 0;
    block_exists_ = false;
  }

#else

  for (size_t i=0; i<block_.size(); i++) {
    delete block_[i];
    block_[i] = 0;
  }

#endif
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
void Patch::p_refresh ()
{
  block_array_->p_refresh(); 
}
#endif

//----------------------------------------------------------------------
#ifdef CONFIG_USE_CHARM
void Patch::p_compute (int cycle, double time, double dt)
{
  DEBUG3 ("cycle %d time %f dt %f",cycle,time,dt);

  block_array_->p_compute(cycle,time,dt);
}
#endif

//========================================================================

#ifndef CONFIG_USE_CHARM
size_t Patch::num_local_blocks() const  throw()
{
  int rank = group_process_->rank();
  return layout_.local_count(rank);
}
#endif

//----------------------------------------------------------------------

#ifndef CONFIG_USE_CHARM
Block * Patch::local_block(size_t i) const throw()
{
  return (i < block_.size()) ? block_[i] : 0;

}
#endif


//======================================================================

void Patch::allocate_array_
(
 bool allocate_blocks,
 const FieldDescr * field_descr
) throw()
  // NOTE: field_descr only needed for MPI; may be null for CHARM++
{

  DEBUG2("patch %p ID = %d",this,id_);
#ifndef CONFIG_USE_CHARM

  // determine local block count nb
  int nb = num_local_blocks();

  // create local blocks
  block_.resize(nb);

#endif

  // Get number of blocks in the patch

  int nbx,nby,nbz;
  layout_.block_count (&nbx, &nby, &nbz);

  // determine block size
  int mbx = size_[0] / nbx;
  int mby = size_[1] / nby;
  int mbz = size_[2] / nbz;

  // Check that blocks evenly subdivide patch
  if (! ((nbx*mbx == size_[0]) &&
	 (nby*mby == size_[1]) &&
	 (nbz*mbz == size_[2]))) {

    ERROR6("Patch::allocate_array",  
	   "Blocks must evenly subdivide Patch: "
	   "patch size = (%d %d %d)  block count = (%d %d %d)",
	   size_[0],size_[1],size_[2],
	   nbx,nby,nbz);
      
  }

  // Determine size of each block
  double xb = (upper_[0] - lower_[0]) / nbx;
  double yb = (upper_[1] - lower_[1]) / nby;
  double zb = (upper_[2] - lower_[2]) / nbz;

  // CREATE AND INITIALIZE NEW DATA BLOCKS

  int num_field_blocks = 1;

#ifdef CONFIG_USE_CHARM

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  DEBUG1 ("simulation = %p",simulation);

  const Factory * factory = simulation->factory();

  block_array_ = new CProxy_Block;

  DEBUG1("ID = %d",id_);
  (*block_array_) = factory->create_block_array
    (nbx,nby,nbz,
     mbx,mby,mbz,
     lower_[0],lower_[1],lower_[2],
     xb,yb,zb,
     thisProxy,
     id_,
     CkMyPe(),
     num_field_blocks,
     allocate_blocks);
    
  block_exists_ = allocate_blocks;
  block_counter_.set_max(nbx*nby*nbz);


#else

  int ip = group_process_->rank();

  for (int ib=0; ib<nb; ib++) {

    // Get index of block ib in the patch

    int ibx,iby,ibz;
    layout_.block_indices (ip,ib, &ibx, &iby, &ibz);

    // create a new data block

    DEBUG1("Patch::Patch(%d)",id_);
    Block * block = factory_->create_block 
      (ibx,iby,ibz,
       nbx,nby,nbz,
       mbx,mby,mbz,
       lower_[0],lower_[1],lower_[2],
       xb,yb,zb,
       id_,
       ip,
       num_field_blocks);

    // Store the data block in the block array
    block_[ib] = block;

    // Allocate data on the block

    block->allocate(field_descr);

  }
#endif /* ! CONFIG_USE_CHARM */
}

//----------------------------------------------------------------------
#include "test.hpp"
 
// Set Patch size, offset, and blocking

static const int patch_size[] = {12,12,12};

static const int patch_offset[] = {5, 2, 9};

static const int patch_blocking[] = {3,3,3};

// Set domain extents

static const double domain_lower[] = {0.0, 0.0, 0.0};
static const double domain_upper[] = {1.0, 1.0, 1.0};

void Patch::p_test () throw()
{

#ifdef CONFIG_USE_CHARM
  Patch * patch = thisProxy.ckLocal();
#else
  Patch * patch = this;
#endif


  unit_assert(patch != NULL);

  //--------------------------------------------------

  unit_func("size");

  int npx,npy,npz;

  patch->size(&npx,&npy,&npz);

  unit_assert(patch_size[0]==npx && 
	      patch_size[1]==npy && 
	      patch_size[2]==npz);

  //--------------------------------------------------

  unit_func("offset");

  int nx0,ny0,nz0;

  patch->offset(&nx0,&ny0,&nz0);

  unit_assert(patch_offset[0]==nx0 && 
	      patch_offset[1]==ny0 && 
	      patch_offset[2]==nz0);

  //--------------------------------------------------

  unit_func("blocking");

  int nbx,nby,nbz;

  patch->blocking(&nbx,&nby,&nbz);

  unit_assert(patch_blocking[0]==nbx && 
	      patch_blocking[1]==nby && 
	      patch_blocking[2]==nbz);

  //--------------------------------------------------

  unit_func("set_lower");

  double xm,ym,zm;
  patch->lower(&xm,&ym,&zm);

  unit_assert(xm==domain_lower[0]);
  unit_assert(ym==domain_lower[1]);
  unit_assert(zm==domain_lower[2]);

  //--------------------------------------------------

  unit_func("set_upper");

  double xp,yp,zp;
  patch->upper(&xp,&yp,&zp);

  unit_assert(xp==domain_upper[0]);
  unit_assert(yp==domain_upper[1]);
  unit_assert(zp==domain_upper[2]);

  //--------------------------------------------------

  unit_func("layout");

  Layout * layout = patch->layout();

  unit_assert(layout != NULL);

  layout->set_process_range(0,1);

  //--------------------------------------------------

#ifdef CONFIG_USE_CHARM
  unit_func("num_blocks");
  unit_assert(patch->num_blocks() == (size_t)nbx*nby*nbz);
#else
  unit_func("num_local_blocks");
  unit_assert(patch->num_local_blocks()==(size_t)nbx*nby*nbz);
#endif

  ItBlock itBlock (patch);

  Block *  block = 0;
  FieldBlock * field_block = 0;

  size_t block_counter = 0;

  const GroupProcess * group_process = GroupProcess::create();

  while ((block = ++itBlock)) {

  //--------------------------------------------------

    unit_func("allocate_blocks");
    unit_assert_quiet(block != NULL);

  //--------------------------------------------------

    unit_class("Block");
    unit_func("field_block");

    field_block = block ? block->field_block() : NULL;

    unit_assert_quiet(field_block != NULL);

    if (block && field_block) {

      //--------------------------------------------------

      unit_class("FieldBlock");
      unit_func("size");

      int nfx, nfy, nfz;

      field_block->size(&nfx,&nfy,&nfz);

      unit_assert_quiet (nfx == patch_size[0] / patch_blocking[0]);
      unit_assert_quiet (nfy == patch_size[1] / patch_blocking[1]);
      unit_assert_quiet (nfz == patch_size[2] / patch_blocking[2]);

      Layout      * layout = patch->layout();

      int ip = group_process->rank();

      int index_local = block_counter;
      int index_global = layout->global_index(ip,index_local);

      //--------------------------------------------------

      unit_class("Layout");
      unit_func("block_indices");

      int ibx,iby,ibz;
      layout->block_indices(index_global,&ibx,&iby,&ibz);

      // Not terribly rigorous
      unit_assert (0 <= ibx && ibx < nbx);
      unit_assert (0 <= iby && iby < nby);
      unit_assert (0 <= ibz && ibz < nbz);

      //--------------------------------------------------

      double xmb,ymb,zmb;
      block->lower (&xmb,&ymb,&zmb);
      double xpb,ypb,zpb;
      block->upper (&xpb,&ypb,&zpb);

      unit_class("Block");
      unit_func("lower");

      unit_assert(cello::err_abs(xm + ibx*(xpb-xmb) , xmb) < 1e-6);
      unit_assert(cello::err_abs(ym + iby*(ypb-ymb) , ymb) < 1e-6);
      unit_assert(cello::err_abs(zm + ibz*(zpb-zmb) , zmb) < 1e-6);

      unit_class("Block");
      unit_func("upper");

      unit_assert(cello::err_abs(xm + (ibx+1)*(xpb-xmb) , xpb) < 1e-6);
      unit_assert(cello::err_abs(ym + (iby+1)*(ypb-ymb) , ypb) < 1e-6);
      unit_assert(cello::err_abs(zm + (ibz+1)*(zpb-zmb) , zpb) < 1e-6);

      //--------------------------------------------------

    }

    //--------------------------------------------------

    ++block_counter;

  }

  //--------------------------------------------------

  unit_class("Block");
  unit_func("index_patch");

  int * b = new int [nbx*nby*nbz];
  int i;
  for (i=0; i<nbx*nby*nbz; i++) b[i]=0;

  while ((block = ++itBlock)) {
    int ibx,iby,ibz;
    block->index_patch(&ibx,&iby,&ibz);
    b[ibx + nbx*(iby + nby*ibz)] = 1;
  }
  for (i=0; i<nbx*nby*nbz; i++) {
    unit_assert_quiet(b[i]==1);
  }

  delete [] b;

  //--------------------------------------------------

#ifndef CONFIG_USE_CHARM
  unit_func("num_local_blocks");
  unit_assert(block_counter == patch->num_local_blocks());
#endif

  
  //--------------------------------------------------

  unit_func("deallocate_blocks");

  patch->deallocate_blocks();

  //--------------------------------------------------

  unit_finalize();

#ifdef CONFIG_USE_CHARM
  proxy_main.p_exit(CkNumPes());
#endif

  PARALLEL_EXIT;
}

