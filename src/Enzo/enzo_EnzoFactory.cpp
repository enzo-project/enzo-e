// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Enzo] Declaration of the EnzoFactory class

#include "enzo.hpp"

#include "charm_enzo.hpp"

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void EnzoFactory::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Factory::pup(p);
}

#endif

//----------------------------------------------------------------------

IoBlock * EnzoFactory::create_io_block () const throw()
{
  return new IoEnzoBlock;
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

CProxy_CommBlock EnzoFactory::create_block_array
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks,
 bool allocate,
 bool testing
 ) const throw()
{
  TRACE ("EnzoFactory::create_block_array");
  CProxy_EnzoBlock enzo_block_array;

  if (allocate) {

    CProxy_ArrayMap array_map  = CProxy_ArrayMap::ckNew(nbx,nby,nbz);
    CkArrayOptions opts;
    opts.setMap(array_map);
    enzo_block_array = CProxy_CommBlock::ckNew(opts);

    //    enzo_block_array = CProxy_EnzoBlock::ckNew();

    int level;

    for (int ix=0; ix<nbx; ix++) {
      for (int iy=0; iy<nby; iy++) {
	for (int iz=0; iz<nbz; iz++) {

	  Index index(ix,iy,iz);

	  TRACE3 ("inserting %d %d %d",ix,iy,iz);
	  enzo_block_array[index].insert 
	    (index,
	     nx,ny,nz,
	     level=0,
	     num_field_blocks,
	     0);

	}
      }
    }

    enzo_block_array.doneInserting();

  } else {

    enzo_block_array = CProxy_EnzoBlock::ckNew();

  }
  TRACE1("EnzoFactory::create_block_array = %p",&enzo_block_array);
  return enzo_block_array;
}

#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

CommBlock * EnzoFactory::create_block
(
#ifdef CONFIG_USE_CHARM
 CProxy_CommBlock block_array,
#else
 Simulation * simulation,
#endif /* CONFIG_USE_CHARM */
 Index index,
 int nx, int ny, int nz,
 int level,
 int num_field_blocks,
 int count_adapt,
 bool testing
 ) const throw()
{

  TRACE("new EnzoBlock");

#ifdef CONFIG_USE_CHARM

   block_array[index].insert
     (
      index,
      nx,ny,nz,
      level,
      num_field_blocks,
      count_adapt,
      testing);

   CommBlock * block = block_array[index].ckLocal();
   TRACE1("block = %p",block);
   //  ASSERT("Factory::create_block()","block is NULL",block != NULL);

   return block;

#else /* CONFIG_USE_CHARM */

  EnzoBlock * enzo_block = new EnzoBlock 
    (
     simulation,
     index,
     nx,ny,nz,
     level,
     num_field_blocks,
     count_adapt);

  return enzo_block;
  
#endif /* CONFIG_USE_CHARM */

}

