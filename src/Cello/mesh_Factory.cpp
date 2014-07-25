// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

#include "mesh.hpp"

//----------------------------------------------------------------------

Hierarchy * Factory::create_hierarchy 
(
 int rank, int refinement,
 int process_first, int process_last_plus) const throw ()
{
  return new Hierarchy 
    (this,rank,refinement,process_first, process_last_plus); 
}

//----------------------------------------------------------------------

void Factory::pup (PUP::er &p)

{
  TRACEPUP;

  PUP::able::pup(p);

  // NOTE: change this function whenever attributes change
}

//----------------------------------------------------------------------

IoBlock * Factory::create_io_block () const throw()
{
  return new IoBlock;
}

//----------------------------------------------------------------------

IoFieldBlock * Factory::create_io_field_block () const throw()
{
  return new IoFieldBlock;

}

//----------------------------------------------------------------------

CProxy_CommBlock Factory::create_block_array
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks,
 bool testing
 ) const throw()
{
  TRACE7("Factory::create_block_array(na(%d %d %d) n(%d %d %d num_field_blocks %d",
	 nbx,nby,nbz,nx,ny,nz,num_field_blocks);

  CProxy_CommBlock proxy_block;

  // --------------------------------------------------
  // ENTRY: #1 Factory::create_block_array() -> ArrayMap::ArrayMap()
  // ENTRY: create
  // --------------------------------------------------
  CProxy_ArrayMap array_map  = CProxy_ArrayMap::ckNew(nbx,nby,nbz);
  // --------------------------------------------------

  CkArrayOptions opts;
  opts.setMap(array_map);
  proxy_block = CProxy_CommBlock::ckNew(opts);

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

	// --------------------------------------------------
	// ENTRY: #2 Factory::create_block_array() -> CommBlock::CommBlock()
	// ENTRY: level == 0 block array insert
	// --------------------------------------------------
	proxy_block[index].insert
	  (index,
	   nx,ny,nz,
	   num_field_blocks,
	   count_adapt = 0,
	   cycle,time,dt,
	   0,NULL,op_array_copy,
	   num_face_level, face_level,
	   testing);
	// --------------------------------------------------

      }
    }
  }

  proxy_block.doneInserting();


  TRACE1("Factory::create_block_array = %p",&proxy_block);
  return proxy_block;
}

//----------------------------------------------------------------------
CommBlock * Factory::create_block
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

  TRACE3("Factory::create_block(%d %d %d)",nx,ny,nz);
  TRACE2("Factory::create_block(num_field_blocks %d  count_adapt %d)",
	 num_field_blocks,count_adapt);

  // --------------------------------------------------
  // ENTRY: #3 Factory::create_block() -> CommBlock::CommBlock()
  // ENTRY: level > 0 block array insert
  // --------------------------------------------------
  (*block_array)[index].insert
    (
     index,
     nx,ny,nz,
     num_field_blocks,
     count_adapt,
     cycle, time,dt,
     narray, array,op_array,
     num_face_level, face_level,
     testing);
  // --------------------------------------------------

  CommBlock * block = (*block_array)[index].ckLocal();

  ASSERT("Factory::create_block()","block is NULL",block != NULL);

  return block;

}

