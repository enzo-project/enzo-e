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
 double xm, double ym, double zm,
 double hx, double hy, double hz,
 int num_field_blocks,
 bool allocate,
 bool testing
 ) const throw()
{
  TRACE ("EnzoFactor::create_block_array");
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
	  index.set_array(ix,iy,iz);
	  index.set_level(0);

	  TRACE3 ("inserting %d %d %d",ix,iy,iz);
	  enzo_block_array[index].insert 
	    (ix,iy,iz,
	     nbx,nby,nbz,
	     nx,ny,nz,
	     level=0,
	     xm,ym,zm, 
	     hx,hy,hz, 
	     num_field_blocks);

	}
      }
    }

    TRACE ("done inserting");
    enzo_block_array.doneInserting();

  } else {

    enzo_block_array = CProxy_EnzoBlock::ckNew();

  }
  TRACE1("EnzoFactory::create_block_array = %p",&enzo_block_array);
  return enzo_block_array;
}

#else

//----------------------------------------------------------------------

CommBlock * EnzoFactory::create_block
(
 int ibx, int iby, int ibz,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double hx, double hy, double hz,
 int num_field_blocks,
 bool testing
 ) const throw()
{
  int level;
#ifdef CONFIG_USE_CHARM
    CProxy_CommBlock block_array = CProxy_EnzoBlock::ckNew
    (nbx,nby,nbz,
     nx,ny,nz,
     level=0,
     xm,ym,zm, 
     xb,yb,zb, 
     num_field_blocks,
     nbx,nby,nbz);

    Index index(ibx,iby,ibz);

    return block_array[index].ckLocal();
   
#else /* CONFIG_USE_CHARM */
  return new EnzoBlock 
    (
     ibx,iby,ibz, 
     nbx,nby,nbz,
     nx,ny,nz,
     level=0,
     xm,ym,zm, 
     hx,hy,hz, 
     num_field_blocks);
#endif /* CONFIG_USE_CHARM */
}

#endif


