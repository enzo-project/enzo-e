// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

#include "mesh.hpp"

//----------------------------------------------------------------------

Hierarchy * Factory::create_hierarchy 
(
#ifndef CONFIG_USE_CHARM
 Simulation * simulation,
#endif
 int dimension, int refinement,
 int process_first, int process_last_plus) const throw ()
{
  return new Hierarchy 
    (
#ifndef CONFIG_USE_CHARM
     simulation,
#endif
     this,dimension,refinement,process_first, process_last_plus); 
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Factory::pup (PUP::er &p)

{
  TRACEPUP;

  PUP::able::pup(p);

  // NOTE: change this function whenever attributes change
}

#endif

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
#ifdef CONFIG_USE_CHARM

CProxy_CommBlock Factory::create_block_array
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks,
 bool allocate,
 bool testing
 ) const throw()
{
  TRACE8("Factory::create_block_array(na(%d %d %d) n(%d %d %d num_field_blocks %d  allocate %d",
	 nbx,nby,nbz,nx,ny,nz,num_field_blocks,allocate);

  CProxy_CommBlock proxy_block;

  if (allocate) {

    CProxy_ArrayMap array_map  = CProxy_ArrayMap::ckNew(nbx,nby,nbz);
    CkArrayOptions opts;
    opts.setMap(array_map);
    proxy_block = CProxy_CommBlock::ckNew(opts);

    int count_adapt;
    bool initial;

    for (int ix=0; ix<nbx; ix++) {
      for (int iy=0; iy<nby; iy++) {
	for (int iz=0; iz<nbz; iz++) {

	  Index index(ix,iy,iz);

	  proxy_block[index].insert 
	    (index,
	     nx,ny,nz,
	     num_field_blocks,
	     count_adapt = 0,
	     initial = true,
	     0,NULL,
	     testing);

	}
      }
    }

    proxy_block.doneInserting();

  } else {

    proxy_block = CProxy_CommBlock::ckNew();

  }

  TRACE1("Factory::create_block_array = %p",&proxy_block);
  return proxy_block;
}

#endif

//----------------------------------------------------------------------
CommBlock * Factory::create_block
(
#ifdef CONFIG_USE_CHARM
 CProxy_CommBlock * block_array,
#else /* CONFIG_USE_CHARM */
 Simulation * simulation,
#endif /* CONFIG_USE_CHARM */
 Index index,
 int nx, int ny, int nz,
 int num_field_blocks,
 int count_adapt,
 bool initial,
 int narray,
 char * array,
 bool testing
 ) const throw()
{

  TRACE6("Factory::create_block(n(%d %d %d)  num_field_blocks %d  count_adatp %d  initial %d)",
	 nx,ny,nz,num_field_blocks,count_adapt,initial);

#ifdef CONFIG_USE_CHARM

  (*block_array)[index].insert
    (
     index,
     nx,ny,nz,
     num_field_blocks,
     count_adapt,
     initial,
     narray,
     array,
     testing);

  CommBlock * block = (*block_array)[index].ckLocal();

  ASSERT("Factory::create_block()","block is NULL",block != NULL);

  return block;

#else /* CONFIG_USE_CHARM */

   CommBlock * comm_block = new CommBlock 
     (simulation,
      index,
      nx,ny,nz,
      num_field_blocks,
      count_adapt,
      initial,
      narray,
      array,
      testing);

   return comm_block;

#endif /* CONFIG_USE_CHARM */
}

