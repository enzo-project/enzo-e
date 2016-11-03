// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Enzo] Declaration of the EnzoFactory class

#include "enzo.hpp"

#include "charm_enzo.hpp"

// #define DEBUG_ENZO_FACTORY

//----------------------------------------------------------------------

void EnzoFactory::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Factory::pup(p);
}

//----------------------------------------------------------------------

IoBlock * EnzoFactory::create_io_block () const throw()
{
  return new IoEnzoBlock;
}

//----------------------------------------------------------------------

CProxy_Block EnzoFactory::create_block_array
(
 DataMsg * data_msg,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks
 ) const throw()
{
  CProxy_EnzoBlock enzo_block_array;

#ifdef DEBUG_ENZO_FACTORY
  CkPrintf ("%d EnzoFactory::create_block_array %d %d %d)\n",CkMyPe(), nx,ny,nz);
  CkPrintf ("%d EnzoFactory::create_block_array() num_field_blocks %d\n",
	    CkMyPe(), num_field_blocks);
#endif
  
  //  CProxy_MappingTree array_map  = CProxy_MappingTree::ckNew(nbx,nby,nbz);
  CProxy_MappingArray array_map  = CProxy_MappingArray::ckNew(nbx,nby,nbz);


  CkArrayOptions opts;
  opts.setMap(array_map);
  TRACE_CHARM("ckNew(nbx,nby,nbz)");
  enzo_block_array = CProxy_EnzoBlock::ckNew(opts);

  int count_adapt;

  int    cycle = 0;
  double time  = 0.0;
  double dt    = 0.0;
  int num_face_level = 0;
  int * face_level = 0;

  for (int ix=0; ix<nbx; ix++) {
    for (int iy=0; iy<nby; iy++) {
      for (int iz=0; iz<nbz; iz++) {

	Index index(ix,iy,iz);

	TRACE3 ("inserting %d %d %d",ix,iy,iz);

	MsgRefine * msg = new MsgRefine 
	  (index,
	   nx,ny,nz,
	   num_field_blocks,
	   count_adapt = 0,
	   cycle,time,dt,
	   refresh_same,
	   num_face_level, face_level);

	msg->set_data_msg(data_msg);

	enzo_block_array[index].insert (msg);

	// --------------------------------------------------

      }
    }
  }

  TRACE1("EnzoFactory::create_block_array = %p",&enzo_block_array);
  return enzo_block_array;
}

//----------------------------------------------------------------------

void EnzoFactory::create_subblock_array
(
 DataMsg * data_msg,
 CProxy_Block * block_array,
 int min_level,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks) const throw()
{
  TRACE8("EnzoFactory::create_subblock_array(min_level %d na(%d %d %d) n(%d %d %d) num_field_blocks %d",
	 min_level,nbx,nby,nbz,nx,ny,nz,num_field_blocks);

  if (min_level >= 0) {
    WARNING1("EnzoFactor::create_subblock_array",
	     "Trying to create subblock array with min_level %d >= 0",
	     min_level);
  }

  CProxy_EnzoBlock * enzo_block_array = 
    static_cast<CProxy_EnzoBlock*> (block_array);

  for (int level = -1; level >= min_level; level--) {

    if (nbx > 1) nbx = ceil(0.5*nbx);
    if (nby > 1) nby = ceil(0.5*nby);
    if (nbz > 1) nbz = ceil(0.5*nbz);

    // CProxy_MappingTree array_map  = CProxy_MappingTree::ckNew(nbx,nby,nbz);
    CProxy_MappingArray array_map  = CProxy_MappingArray::ckNew(nbx,nby,nbz);

    CkArrayOptions opts;
    opts.setMap(array_map);
    TRACE_CHARM("ckNew(nbx,nby,nbz)");

    int count_adapt;

    int    cycle = 0;
    double time  = 0.0;
    double dt    = 0.0;
    int num_face_level = 0;
    int * face_level = 0;

    for (int ix=0; ix<nbx; ix++) {
      for (int iy=0; iy<nby; iy++) {
	for (int iz=0; iz<nbz; iz++) {

	  Index index(ix,iy,iz);

	  index.set_level(level);

	  TRACE3 ("inserting %d %d %d",ix,iy,iz);

	  MsgRefine * msg = new MsgRefine 
	    (index,
	     nx,ny,nz,
	     num_field_blocks,
	     count_adapt=0,
	     cycle,time,dt,
	     refresh_same,
	     num_face_level, face_level);

	  msg->set_data_msg(data_msg);

	  (*enzo_block_array)[index].insert (msg);

	  // --------------------------------------------------

	}
      }
    }

  }

}

//----------------------------------------------------------------------

Block * EnzoFactory::create_block
(
 DataMsg * data_msg,
 CProxy_Block * block_array,
 Index index,
 int nx, int ny, int nz,
 int num_field_blocks,
 int count_adapt,
 int cycle, double time, double dt,
 int narray, char * array, int refresh_type,
 int num_face_level, int * face_level,
 Simulation * simulation
 ) const throw()
{
#ifdef DEBUG_ENZO_FACTORY
  CkPrintf ("%d EnzoFactory::create_block %d %d %d)\n",CkMyPe(), nx,ny,nz);
  CkPrintf ("%d EnzoFactory::create_block() num_field_blocks %d  count_adapt %d\n",
	    CkMyPe(), num_field_blocks,count_adapt);
#endif

  TRACE3("EnzoFactory::create_block(%d %d %d)",nx,ny,nz);
  TRACE2("EnzoFactory::create_block() num_field_blocks %d  count_adapt %d",
	 num_field_blocks,count_adapt);

  CProxy_EnzoBlock * enzo_block_array = (CProxy_EnzoBlock * ) block_array;

#ifdef CELLO_DEBUG
  index.print("DEBUG insert()",-1,2,false,simulation);
#endif
  
  MsgRefine * msg = new MsgRefine 
    (index,
     nx,ny,nz,
     num_field_blocks,
     count_adapt,
     cycle,time,dt,
     refresh_type,
     num_face_level, face_level);

  msg->set_data_msg(data_msg);

  (*enzo_block_array)[index].insert ( msg );

  // --------------------------------------------------

#ifdef CELLO_TRACE
  index.print("ADAPT REFINE insert()",-1,2,false,simulation);
#endif
  EnzoBlock * block = (*enzo_block_array)[index].ckLocal();
  TRACE1("block = %p",block);
  //  ASSERT("Factory::create_block()","block is NULL",block != NULL);

  return block;

}

