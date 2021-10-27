// See LICENSE_CELLO file for license and copyright information

/// @file     test_Adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-10-27
/// @brief    Test program for the Adapt class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Adapt");

  // 4            ^
  // 3          o o
  // 2    o   o     o
  // 1  o   o         o
  // 0  V   V
  //
  Index index[8] =
    {
     Index(0,0,0),
     Index(0,0,1),
     Index(0,1,0),
     Index(0,1,1),
     Index(1,0,0),
     Index(1,0,1),
     Index(1,1,0),
     Index(1,1,1) };

  const int level_min[8] = {0,2,0,2,3,4,2,1};
  const int level_max[8] = {2,3,2,3,4,4,3,2};
  Adapt adapt[8];

  // Initialize level bounds for each "block"
  for (int i0=0; i0<8; i0++) {
    int ip = (i0 + 1) % 8;
    int im = (i0 + 8 - 1) % 8;
    adapt[i0].allocate_level_bounds(2);
    adapt[i0].set_initial_level_bounds
      (0,index[i0],level_min[i0], level_max[i0]);
    adapt[i0].set_initial_level_bounds
      (1,index[im],level_min[im], level_max[im]);
    adapt[i0].set_initial_level_bounds
      (2,index[ip],level_min[ip], level_max[ip]);
  }

  
  unit_func("is_committed()");

  unit_assert (adapt[0].is_committed() == false);
  unit_assert (adapt[1].is_committed() == false);
  unit_assert (adapt[2].is_committed() == false);
  unit_assert (adapt[3].is_committed() == false);
  unit_assert (adapt[4].is_committed() == false);
  unit_assert (adapt[5].is_committed() == true);
  unit_assert (adapt[6].is_committed() == false);
  unit_assert (adapt[7].is_committed() == false);

  // Update neighbor level bounds

  CkPrintf ("committed: ");
  for (int i0=0; i0<8; i0++) CkPrintf ("%1d",adapt[i0].is_committed());
  CkPrintf ("\n");
  for (int i0=0; i0<8; i0++) {
    int min,max;
    int im = (i0 + 8 - 1) % 8;
    adapt[im].get_level_bounds(&min,&max);
    adapt[i0].update_level_bounds(index[im],min,max);
    int ip = (i0 + 8 + 1) % 8;
    adapt[ip].get_level_bounds(&min,&max);
    adapt[i0].update_level_bounds(index[ip],min,max);
    adapt[i0].evaluate_level_bounds();
  }
  CkPrintf ("committed: ");
  for (int i0=0; i0<8; i0++) CkPrintf ("%1d",adapt[i0].is_committed());
  CkPrintf ("\n");
  for (int i0=0; i0<8; i0++) {
    int min,max;
    int im = (i0 + 8 - 1) % 8;
    adapt[im].get_level_bounds(&min,&max);
    adapt[i0].update_level_bounds(index[im],min,max);
    int ip = (i0 + 8 + 1) % 8;
    adapt[ip].get_level_bounds(&min,&max);
    adapt[i0].update_level_bounds(index[ip],min,max);
    adapt[i0].evaluate_level_bounds();
  }
  CkPrintf ("committed: ");
  for (int i0=0; i0<8; i0++) CkPrintf ("%1d",adapt[i0].is_committed());
  CkPrintf ("\n");
  for (int i0=0; i0<8; i0++) {
    int min,max;
    int im = (i0 + 8 - 1) % 8;
    adapt[im].get_level_bounds(&min,&max);
    adapt[i0].update_level_bounds(index[im],min,max);
    int ip = (i0 + 8 + 1) % 8;
    adapt[ip].get_level_bounds(&min,&max);
    adapt[i0].update_level_bounds(index[ip],min,max);
    adapt[i0].evaluate_level_bounds();
  }
  CkPrintf ("committed: ");
  for (int i0=0; i0<8; i0++) CkPrintf ("%1d",adapt[i0].is_committed());
  CkPrintf ("\n");
  // adapt[0].update_level_bounds (index[2],3,3);
  // adapt[0].update_level_bounds (index[6],3,3);
  // unit_func("is_committed()");
  // unit_assert (adapt[0].is_committed(0) == false);
  // unit_assert (adapt[0].is_committed(1) == false);
  // unit_assert (adapt[0].is_committed(2) == true);
  // unit_assert (adapt[0].is_committed(3) == true);
  // unit_assert (adapt[0].is_committed(4) == true);
  // unit_assert (adapt[0].is_committed(5) == false);
  // unit_assert (adapt[0].is_committed(6) == true);
  // unit_assert (adapt[0].is_committed(7) == false);
  // unit_assert (adapt[0].is_committed(adapt[0].index(index[0])) == false);
  // unit_assert (adapt[0].is_committed(adapt[0].index(index[1])) == false);
  // unit_assert (adapt[0].is_committed(adapt[0].index(index[2])) == true);
  // unit_assert (adapt[0].is_committed(adapt[0].index(index[3])) == true);
  // unit_assert (adapt[0].is_committed(adapt[0].index(index[4])) == true);
  // unit_assert (adapt[0].is_committed(adapt[0].index(index[5])) == false);
  // unit_assert (adapt[0].is_committed(adapt[0].index(index[6])) == true);
  // unit_assert (adapt[0].is_committed(adapt[0].index(index[7])) == false);

  // unit_func("evaluate_level_bounds()");
  
  // unit_assert (adapt[0].evaluate_level_bounds());
  // unit_assert (! adapt[0].evaluate_level_bounds());

  // int min,max;

  // unit_func("get_level_bounds()");

  // adapt[0].get_level_bounds (&min,&max);

  // unit_assert (min==2);
  // unit_assert (max==3);
  //--------------------------------------------------

  unit_func ("function()");

  //--------------------------------------------------

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

