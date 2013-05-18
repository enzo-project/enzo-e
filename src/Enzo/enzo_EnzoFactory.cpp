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
  TRACE8("EnzoFactory::create_block_array(na(%d %d %d) n(%d %d %d num_field_blocks %d  allocate %d",
	 nbx,nby,nbz,nx,ny,nz,num_field_blocks,allocate);

  CProxy_EnzoBlock enzo_block_array;

  if (allocate) {

    CProxy_ArrayMap array_map  = CProxy_ArrayMap::ckNew(nbx,nby,nbz);
    CkArrayOptions opts;
    opts.setMap(array_map);
    enzo_block_array = CProxy_CommBlock::ckNew(opts);

    int level;

    int count_adapt;

    bool initial;

    int    cycle = 0;
    double time  = 0.0;
    double dt    = 0.0;

    for (int ix=0; ix<nbx; ix++) {
      for (int iy=0; iy<nby; iy++) {
	for (int iz=0; iz<nbz; iz++) {

	  Index index(ix,iy,iz);

	  TRACE3 ("inserting %d %d %d",ix,iy,iz);
	  enzo_block_array[index].insert 
	    (index,
	     nx,ny,nz,
	     num_field_blocks,
	     count_adapt = 0,
	     initial=true,
	     cycle, time, dt,
	     0,NULL,op_array_copy,
	     testing);

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
 CProxy_CommBlock * block_array,
#else
 Simulation * simulation,
#endif /* CONFIG_USE_CHARM */
 Index index,
 int nx, int ny, int nz,
 int num_field_blocks,
 int count_adapt,
 bool initial,
 int cycle, double time, double dt,
 int narray, char * array, int op_array,
 bool testing
 ) const throw()
{
  TRACE6("EnzoFactory::create_block(n(%d %d %d)  num_field_blocks %d  count_adatp %d  initial %d)",
	 nx,ny,nz,num_field_blocks,count_adapt,initial);


#ifdef CONFIG_USE_CHARM

  CProxy_EnzoBlock * enzo_block_array = (CProxy_EnzoBlock * ) block_array;

  (*enzo_block_array)[index].insert
     (
      index,
      nx,ny,nz,
      num_field_blocks,
      count_adapt,
      initial,
      cycle,time,dt,
      narray, array, op_array,
      testing);

  CommBlock * block = (*enzo_block_array)[index].ckLocal();
   TRACE1("block = %p",block);
   //  ASSERT("Factory::create_block()","block is NULL",block != NULL);

   return block;

#else /* CONFIG_USE_CHARM */

  EnzoBlock * enzo_block = new EnzoBlock 
    (
     simulation,
     index,
     nx,ny,nz,
     num_field_blocks,
     initial,
     count_adapt,
     cycle, time, dt,
     narray, array, op_array,
     testing);

  return enzo_block;
  
#endif /* CONFIG_USE_CHARM */

}

