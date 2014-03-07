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

CProxy_CommBlock EnzoFactory::create_block_array
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks,
 bool testing
 ) const throw()
{
  TRACE7("EnzoFactory::create_block_array(na(%d %d %d) n(%d %d %d num_field_blocks %d",	 nbx,nby,nbz,nx,ny,nz,num_field_blocks);

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

	// --------------------------------------------------
	// ENTRY: #2 EnzoFactory::create_block_array() -> EnzoBlock::EnzoBlock()
	// ENTRY: level == 0 block array insert
	// --------------------------------------------------
	enzo_block_array[index].insert 
	  (index,
	   nx,ny,nz,
	   num_field_blocks,
	   count_adapt = 0,
	   cycle, time, dt,
	   0,NULL,op_array_copy,
	   num_face_level, face_level,
	   testing);
	// --------------------------------------------------

      }
    }
  }

  enzo_block_array.doneInserting();


  TRACE1("EnzoFactory::create_block_array = %p",&enzo_block_array);
  return enzo_block_array;
}

//----------------------------------------------------------------------

CommBlock * EnzoFactory::create_block
(
 CProxy_CommBlock * block_array,
 Index index,
 int nx, int ny, int nz,
 int num_field_blocks,
 int count_adapt,
 int cycle, double time, double dt,
 int narray, char * array, int op_array,
 int num_face_level, int * face_level,
 bool testing,
 Simulation * simulation
 ) const throw()
{
  TRACE3("EnzoFactory::create_block(%d %d %d)",nx,ny,nz);
  TRACE2("EnzoFactory::create_block(num_field_blocks %d  count_adapt %d",
	 num_field_blocks,count_adapt);


  CProxy_EnzoBlock * enzo_block_array = (CProxy_EnzoBlock * ) block_array;

#ifdef CELLO_DEBUG
  index.print("DEBUG insert()",-1,2,false,simulation);
#endif

  // --------------------------------------------------
  // ENTRY: #3 EnzoFactory::create_block() -> EnzoBlock::EnzoBlock()
  // ENTRY: level > 0 block array insert
  // --------------------------------------------------
  (*enzo_block_array)[index].insert
    (
     index,
     nx,ny,nz,
     num_field_blocks,
     count_adapt,
     cycle,time,dt,
     narray, array, op_array,
     num_face_level, face_level,
     testing);
  // --------------------------------------------------

#ifdef CELLO_TRACE
  index.print("ADAPT REFINE insert()",-1,2,false,simulation);
#endif
  EnzoBlock * block = (*enzo_block_array)[index].ckLocal();
  TRACE1("block = %p",block);
  //  ASSERT("Factory::create_block()","block is NULL",block != NULL);

  return block;

}

