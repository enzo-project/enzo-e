// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-19
/// @brief    Unit tests for the Layout class

#include <math.h>
#include "cello.hpp"

#include "error.hpp"
#include "test.hpp"
#include "parallel.hpp"

#define TOL 2e-16

#include "parallel.def"
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
  int nb,nbx,nby,nbz;

  //--------------------------------------------------
  // offset 0  count 1  blocks (1,1,1)
  //--------------------------------------------------

  unit_func("processors");

  layout.set_processors(0,1);
  layout.processors(&p0,&np);

  unit_assert (p0 == 0);
  unit_assert (np == 1);
  
  unit_func("blocks");

  layout.set_blocks (1,1,1);
  nb = layout.blocks(&nbx,&nby,&nbz);

  unit_assert (nb == nbx*nby*nbz);
  unit_assert (nbx == 1);
  unit_assert (nby == 1);
  unit_assert (nbz == 1);

  unit_func("process");

  unit_assert(layout.process (0,0,0) == 0);

  //--------------------------------------------------
  // offset 0  count 1  blocks (5,3,2)
  //--------------------------------------------------

  
  unit_func("processors");

  layout.set_processors(0,1);
  layout.processors(&p0,&np);

  unit_assert (p0 == 0);
  unit_assert (np == 1);
  
  unit_func("blocks");

  layout.set_blocks (5,3,2);
  nb = layout.blocks(&nbx,&nby,&nbz);

  unit_assert (nb == nbx*nby*nbz);
  unit_assert (nbx == 5);
  unit_assert (nby == 3);
  unit_assert (nbz == 2);

  unit_func("process");

  unit_assert(layout.process (0,0,0) == 0);
  unit_assert(layout.process (5-1,3-1,2-1) == 0);


  //--------------------------------------------------
  // offset 0  count 30  blocks (5,3,2)
  //--------------------------------------------------
  
  unit_func("processors");

  layout.set_processors(0,30);
  layout.processors(&p0,&np);

  unit_assert (p0 == 0);
  unit_assert (np == 30);
  
  unit_func("blocks");

  layout.set_blocks (5,3,2);
  nb = layout.blocks(&nbx,&nby,&nbz);

  unit_assert (nb == nbx*nby*nbz);
  unit_assert (nbx == 5);
  unit_assert (nby == 3);
  unit_assert (nbz == 2);

  unit_func("process");

  unit_assert(layout.process (0,0,0)       == 0);
  unit_assert(layout.process (5-1,3-1,2-1) == 30-1);

  //--------------------------------------------------
  // offset 7  count 30  blocks (5,3,2)
  //--------------------------------------------------
  
  unit_func("processors");

  layout.set_processors(7,30);
  layout.processors(&p0,&np);

  unit_assert (p0 == 7);
  unit_assert (np == 30);
  
  unit_func("blocks");

  layout.set_blocks (5,3,2);
  nb = layout.blocks(&nbx,&nby,&nbz);

  unit_assert (nb == nbx*nby*nbz);
  unit_assert (nbx == 5);
  unit_assert (nby == 3);
  unit_assert (nbz == 2);

  unit_func("process");

  unit_assert(layout.process (0,0,0)       == 7);
  unit_assert(layout.process (5-1,3-1,2-1) == 7+30-1);


  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Layout.def.h)
