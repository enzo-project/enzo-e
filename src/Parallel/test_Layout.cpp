// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-19
/// @brief    Unit tests for the Layout class

#include "test.hpp"

#include "parallel.hpp"

#include PARALLEL_CHARM_INCLUDE(test_Layout.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  unit_init (parallel->rank(), parallel->size());

  unit_class("Layout");

  unit_func("Layout");
  Layout layout;
  unit_assert(true);

  int p0,np;
  int nbx,nby,nbz,nb;
  int ibx,iby,ibz;

  //--------------------------------------------------
  // first 0  count 1  block_count (1,1,1)
  //--------------------------------------------------

  unit_func("process_range");

  layout.set_process_range(0,1);
  layout.process_range(&p0,&np);

  unit_assert (p0 == 0);
  unit_assert (np == 1);
  
  unit_func("block_count");

  layout.set_block_count (1,1,1);
  nb = layout.block_count(&nbx,&nby,&nbz);

  unit_assert (nb == nbx*nby*nbz);
  unit_assert (nbx == 1);
  unit_assert (nby == 1);
  unit_assert (nbz == 1);

  unit_func("process");

  unit_assert(layout.process (0) == 0);

  unit_assert(layout.process (-1) == PROCESS_NULL);
  unit_assert(layout.process (nb) == PROCESS_NULL);

  unit_func("local_count");

  unit_assert(layout.local_count(0)==nb);

  unit_assert(layout.local_count(-1)==0);
  unit_assert(layout.local_count(np)==0);

  unit_func("block_indices");

  ibx = iby = ibz = -1;

  layout.block_indices(0,&ibx,&iby,&ibz);

  unit_assert(ibx == 0);
  unit_assert(iby == 0);
  unit_assert(ibz == 0);

  unit_func("block_index");

  unit_assert(layout.block_index(ibx,iby,ibz) == 0);


  //--------------------------------------------------
  // first 0  count 1  block_count (5,3,2)
  //--------------------------------------------------

  
  unit_func("process_range");

  layout.set_process_range(0,1);
  layout.process_range(&p0,&np);

  unit_assert (p0 == 0);
  unit_assert (np == 1);
  
  unit_func("block_count");

  layout.set_block_count (5,3,2);
  nb = layout.block_count(&nbx,&nby,&nbz);

  unit_assert (nb == nbx*nby*nbz);
  unit_assert (nbx == 5);
  unit_assert (nby == 3);
  unit_assert (nbz == 2);

  unit_func("process");

  unit_assert(layout.process (layout.block_index(0,0,0)) == 0);
  unit_assert(layout.process (layout.block_index(5-1,3-1,2-1)) == 0);

  unit_assert(layout.process (-1) == PROCESS_NULL);
  unit_assert(layout.process (nb) == PROCESS_NULL);

  unit_func("local_count");

  unit_assert(layout.local_count(0)==nb);

  unit_assert(layout.local_count(-1)==0);
  unit_assert(layout.local_count(np)==0);

  unit_func("block_indices");

  ibx = iby = ibz = -1;

  layout.block_indices(0,&ibx,&iby,&ibz);

  unit_assert(ibx == 0);
  unit_assert(iby == 0);
  unit_assert(ibz == 0);

  unit_func("block_index");

  unit_assert(layout.block_index(ibx,iby,ibz) == 0);

  unit_func("block_indices");

  ibx = iby = ibz = -1;

  layout.block_indices(nb-1,&ibx,&iby,&ibz);

  unit_assert(ibx == 5-1);
  unit_assert(iby == 3-1);
  unit_assert(ibz == 2-1);

  unit_func("block_index");

  unit_assert(layout.block_index(ibx,iby,ibz) == nb-1);


  //--------------------------------------------------
  // first 0  count 30  block_count (5,3,2)
  //--------------------------------------------------
  
  unit_func("process_range");

  layout.set_process_range(0,30);
  layout.process_range(&p0,&np);

  unit_assert (p0 == 0);
  unit_assert (np == 30);
  
  unit_func("block_count");

  layout.set_block_count (5,3,2);
  nb = layout.block_count(&nbx,&nby,&nbz);

  unit_assert (nb == nbx*nby*nbz);
  unit_assert (nbx == 5);
  unit_assert (nby == 3);
  unit_assert (nbz == 2);

  unit_func("process");

  unit_assert(layout.process (layout.block_index(0,0,0))       == 0);
  unit_assert(layout.process (layout.block_index(5-1,3-1,2-1)) == 30-1);

  unit_assert(layout.process (-1) == PROCESS_NULL);
  unit_assert(layout.process (nb) == PROCESS_NULL);

  unit_func("local_count");

  unit_assert(layout.local_count(0)==1);
  unit_assert(layout.local_count(np-1)==1);

  unit_assert(layout.local_count(-1)==0);
  unit_assert(layout.local_count(np)==0);

  unit_func("block_indices");

  ibx = iby = ibz = -1;

  layout.block_indices(0,&ibx,&iby,&ibz);

  unit_assert(ibx == 0);
  unit_assert(iby == 0);
  unit_assert(ibz == 0);

  unit_func("block_index");

  unit_assert(layout.block_index(ibx,iby,ibz) == 0);

  ibx = iby = ibz = -1;

  layout.block_indices(nb-1,&ibx,&iby,&ibz);

  unit_assert(ibx == 5-1);
  unit_assert(iby == 3-1);
  unit_assert(ibz == 2-1);

  unit_func("block_index");

  unit_assert(layout.block_index(ibx,iby,ibz) == nb-1);

  //--------------------------------------------------
  // first 7  count 30  block_count (5,3,2)
  //--------------------------------------------------
  
  unit_func("process_range");

  layout.set_process_range(7,30);
  layout.process_range(&p0,&np);

  unit_assert (p0 == 7);
  unit_assert (np == 30);
  
  unit_func("block_count");

  layout.set_block_count (5,3,2);
  nb = layout.block_count(&nbx,&nby,&nbz);

  unit_assert (nb == nbx*nby*nbz);
  unit_assert (nbx == 5);
  unit_assert (nby == 3);
  unit_assert (nbz == 2);

  unit_func("process");

  unit_assert(layout.process (layout.block_index(0,0,0))       == 7);
  unit_assert(layout.process (layout.block_index(5-1,3-1,2-1)) == 7+30-1);

  unit_assert(layout.process (-1) == PROCESS_NULL);
  unit_assert(layout.process (nb) == PROCESS_NULL);

  unit_func("local_count");

  unit_assert(layout.local_count(p0)==1);
  unit_assert(layout.local_count(p0+np-1)==1);

  unit_assert(layout.local_count(p0-1)==0);
  unit_assert(layout.local_count(p0+np)==0);

  unit_func("block_indices");

  ibx = iby = ibz = -1;

  layout.block_indices(0,&ibx,&iby,&ibz);

  unit_assert(ibx == 0);
  unit_assert(iby == 0);
  unit_assert(ibz == 0);

  unit_func("block_index");

  unit_assert(layout.block_index(ibx,iby,ibz) == 0);

  ibx = iby = ibz = -1;

  layout.block_indices(nb-1,&ibx,&iby,&ibz);

  unit_assert(ibx == 5-1);
  unit_assert(iby == 3-1);
  unit_assert(ibz == 2-1);

  unit_func("block_index");

  unit_assert(layout.block_index(ibx,iby,ibz) == nb-1);

  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Layout.def.h)
