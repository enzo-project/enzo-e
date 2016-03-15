// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Enzo] Declaration of the EnzoFactory class

#include "enzo.hpp"

#include "charm_enzo.hpp"

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
 DataMsgRefine * data_msg,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks,
 bool testing
 ) const throw()
{
  TRACE7("EnzoFactory::create_block_array(na(%d %d %d) n(%d %d %d) num_field_blocks %d",	 nbx,nby,nbz,nx,ny,nz,num_field_blocks);

  CProxy_EnzoBlock enzo_block_array;

  // --------------------------------------------------
  // ENTRY: #1 Factory::create_block_array() -> ArrayMap::ArrayMap()
  // ENTRY: create
  // --------------------------------------------------
  CProxy_ArrayMap array_map  = CProxy_ArrayMap::ckNew(nbx,nby,nbz);
  // --------------------------------------------------

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

#ifdef NEW_REFRESH_REFINE
	enzo_block_array[index].insert (data_msg);
#else
	enzo_block_array[index].insert 
	  (index,
	   nx,ny,nz,
	   num_field_blocks,
	   count_adapt = 0,
	   cycle, time, dt,
	   0,NULL,refresh_same,
	   num_face_level, face_level,
	   testing);
#endif

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
 DataMsgRefine * data_msg,
 CProxy_Block * block_array,
 int min_level,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks,
 bool testing
 ) const throw()
{
  TRACE8("EnzoFactory::create_subblock_array(min_level %d na(%d %d %d) n(%d %d %d) num_field_blocks %d",
	 min_level,nbx,nby,nbz,nx,ny,nz,num_field_blocks);

  if (min_level >= 0) {
    WARNING1("EnzoFactor::create_subblock_array",
	     "Trying to create subblock array with min_level %d >= 0",
	     min_level);
  }

  CProxy_EnzoBlock * enzo_block_array = static_cast<CProxy_EnzoBlock*> (block_array);

  for (int level = -1; level >= min_level; level--) {

    if (nbx > 1) nbx = ceil(0.5*nbx);
    if (nby > 1) nby = ceil(0.5*nby);
    if (nbz > 1) nbz = ceil(0.5*nbz);

    // --------------------------------------------------
    CProxy_ArrayMap array_map  = CProxy_ArrayMap::ckNew(nbx,nby,nbz);
    // --------------------------------------------------

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

#ifdef NEW_REFRESH_REFINE
	  (*enzo_block_array)[index].insert (data_msg);
#else
	  (*enzo_block_array)[index].insert 
	  (index,
	   nx,ny,nz,
	   num_field_blocks,
	   count_adapt = 0,
	   cycle, time, dt,
	   0,NULL,refresh_same,
	   num_face_level, face_level,
	   testing);
#endif
	  // --------------------------------------------------

	}
      }
    }

  }

}

//----------------------------------------------------------------------

Block * EnzoFactory::create_block
(
 DataMsgRefine * data_msg,
 CProxy_Block * block_array,
 Index index,
 int nx, int ny, int nz,
 int num_field_blocks,
 int count_adapt,
 int cycle, double time, double dt,
 int narray, char * array, int refresh_type,
 int num_face_level, int * face_level,
 bool testing,
 Simulation * simulation
 ) const throw()
{
  TRACE3("EnzoFactory::create_block(%d %d %d)",nx,ny,nz);
  TRACE2("EnzoFactory::create_block() num_field_blocks %d  count_adapt %d",
	 num_field_blocks,count_adapt);

  CProxy_EnzoBlock * enzo_block_array = (CProxy_EnzoBlock * ) block_array;

#ifdef CELLO_DEBUG
  index.print("DEBUG insert()",-1,2,false,simulation);
#endif

#ifdef NEW_REFRESH_REFINE
  (*enzo_block_array)[index].insert ( data_msg )
#else
  (*enzo_block_array)[index].insert 
    (     index,
	  nx,ny,nz,
	  num_field_blocks,
	  count_adapt,
	  cycle,time,dt,
	  narray, array, refresh_type,
	  num_face_level, face_level,
	  testing);
#endif
  // --------------------------------------------------

#ifdef CELLO_TRACE
  index.print("ADAPT REFINE insert()",-1,2,false,simulation);
#endif
  EnzoBlock * block = (*enzo_block_array)[index].ckLocal();
  TRACE1("block = %p",block);
  //  ASSERT("Factory::create_block()","block is NULL",block != NULL);

  return block;

}

