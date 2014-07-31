// See LICENSE_CELLO file for license and copyright information

/// @file     test_ItChild.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the ItChild class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

void clear (bool * c, int n) {for (int i=0; i<n; i++) c[i]=0; }

#define IC(ix,iy,iz) ( ((ix+2)%2) + 2*( ((iy+2)%2) + 2*( ((iz+2)%2) )))

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ItChild");

  bool ic[8];
 
  int count;
  int ic3[3];

  //--------------------------------------------------

  // 1D 0
  count = 0;
  unit_func("ItChild(1)");
  clear(ic,2);
  ItChild it_child1 (1); 
  while (it_child1.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 2);
  unit_assert(ic[IC(0,0,0)] && ic[IC(1,0,0)]);

  // 1D -1
  count = 0;
  int if3[3];
  if3[0] = -1;
  ItChild it_child1_m(1,if3);
  clear(ic,2);
  while (it_child1_m.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 1);
  unit_assert(ic[IC(0,0,0)] && ! ic[IC(1,0,0)]);

  // 1D +1
  count = 0;
  if3[0] = 1;
  ItChild it_child1_p(1,if3);
  clear(ic,2);
  while (it_child1_p.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 1);
  unit_assert(! ic[IC(0,0,0)] && ic[IC(1,0,0)]);

  //--------------------------------------------------

  // 2D 0 0

  count = 0;
  unit_func("ItChild(2)");
  clear(ic,4);
  if3[0] = 0;
  if3[1] = 0;
  ItChild it_child2 (2); 
  while (it_child2.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 4);
  unit_assert(ic[IC(0,0,0)] && ic[IC(1,0,0)] && ic[IC(0,1,0)] && ic[IC(1,1,0)]);

  // 2D 0 1

  count = 0;
  unit_func("ItChild(2)");
  clear(ic,4);
  if3[0] = 0;
  if3[1] = 1;
  ItChild it_child2_0p (2,if3); 
  while (it_child2_0p.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 2);
  unit_assert(ic[IC(0,1,0)] && ic[IC(1,1,0)]);


  // 2D 1 -1

  count = 0;
  unit_func("ItChild(2)");
  clear(ic,4);
  if3[0] = 1;
  if3[1] = -1;
  ItChild it_child2_pm (2,if3); 
  while (it_child2_pm.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 1);
  unit_assert(ic[IC(1,0,0)]);


  //--------------------------------------------------

  count = 0;
  unit_func("ItChild(3)");
  clear(ic,8);
  ItChild it_child3 (3); 
  while (it_child3.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 8);
  unit_assert(ic[IC(0,0,0)] && ic[IC(1,0,0)] && ic[IC(0,1,0)] && ic[IC(1,1,0)]
           && ic[IC(0,0,1)] && ic[IC(1,0,1)] && ic[IC(1,1,0)] && ic[IC(1,1,1)]);

  // 2D 0 0 -1

  count = 0;
  unit_func("ItChild(2)");
  clear(ic,4);
  if3[0] = 0;
  if3[1] = 0;
  if3[2] = -1;
  ItChild it_child3_00m (2,if3); 
  while (it_child3_00m.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 4);
  unit_assert(ic[IC(0,0,0)] && ic[IC(1,0,0)] && ic[IC(0,1,0)] && ic[IC(1,1,0)]);


  // 2D 1 0 1

  count = 0;
  unit_func("ItChild(2)");
  clear(ic,4);
  if3[0] = 1;
  if3[1] = 0;
  if3[2] = 1;
  ItChild it_child2_p0p (2,if3); 
  while (it_child2_p0p.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 2);
  unit_assert(ic[5] && ic[7]);


  // 2D -1 -1 -1

  count = 0;
  unit_func("ItChild(2)");
  clear(ic,4);
  if3[0] = -1;
  if3[1] = -1;
  if3[2] = -1;
  ItChild it_child3_mmm (2,if3); 
  while (it_child3_mmm.next(ic3)) {
    ic[IC3(ic3)] = true;
    ++count;
  }
  unit_assert(count == 1);
  unit_assert(ic[IC(0,0,0)]);


  //--------------------------------------------------

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

