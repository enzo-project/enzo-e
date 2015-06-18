// See LICENSE_CELLO file for license and copyright information

/// @file     test_ItFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the ItFace class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

void clear (bool * f, int n) {for (int i=0; i<n; i++) f[i]=0; }

#define IF(ix,iy,iz) ((ix+1) + 3*((iy+1) + 3*(iz+1)))

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ItFace");

  int count;
  int if3[3];
  int ic3[3];
  int ipf3[3];
  bool iface[27];
  
  // -------- 1D --------

  bool periodic[3][2] = {{true,true},{true,true},{true,true}};
  int n3[3] = { 4, 4, 4 };
  Index index;
  count = 0;
  clear(iface,27);
  unit_func("ItFace(1,0)");
  ItFace it_face10 (1,0,periodic,n3,index); 
  while (it_face10.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 3 - 1);
  unit_assert(iface[IF(-1,0,0)] && iface[IF(+1,0,0)]);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(1,0,child(0))");
  ic3[0] = 0;
  ic3[1] = 0;
  ic3[2] = 0;
  ItFace it_face10_0 (1,0,periodic,n3,index,ic3); 
  while (it_face10_0.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 1);
  unit_assert(iface[IF(-1,0,0)]);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(1,0,child(1))");
  ic3[0] = 1;
  ic3[1] = 0;
  ic3[2] = 0;
  ItFace it_face10_1 (1,0,periodic,n3,index,ic3); 
  while (it_face10_1.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 1);
  unit_assert(iface[IF(+1,0,0)]);

  // -------- 2D --------

  count = 0;
  clear(iface,27);
  unit_func("ItFace(2,0)");
  ItFace it_face20 (2,0,periodic,n3,index);
  while (it_face20.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 8);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(2,0,child(0,0))");
  ic3[0] = 0;
  ic3[1] = 0;
  ic3[2] = 0;
  ItFace it_face20_00 (2,0,periodic,n3,index,ic3); 
  while (it_face20_00.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 3);
  unit_assert(iface[IF(-1, 0,0)] && 
	      iface[IF(-1,-1,0)] && 
	      iface[IF( 0,-1,0)]);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(2,0,child(0,1))");
  ic3[0] = 0;
  ic3[1] = 1;
  ic3[2] = 0;
  ItFace it_face20_01 (2,0,periodic,n3,index,ic3); 
  while (it_face20_01.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 3);
  unit_assert(iface[IF(-1, 0, 0)] && 
	      iface[IF(-1, 1, 0)] && 
	      iface[IF( 0, 1, 0)]);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(2,1)");
  ItFace it_face21 (2,1,periodic,n3,index); 
  while (it_face21.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 4);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(2,1,child(1,1))");
  ic3[0] = 1;
  ic3[1] = 1;
  ic3[2] = 0;
  ItFace it_face21_11 (2,1,periodic,n3,index,ic3); 
  while (it_face21_11.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 2);
  unit_assert(iface[IF(0, 1, 0)] && 
	      iface[IF(1, 0, 0)]);


  count = 0;
  clear(iface,27);
  unit_func("ItFace(2,1,child(1,1),face(1,0))");
  ic3[0] = 1;
  ic3[1] = 1;
  ic3[2] = 0;
  ipf3[0] = 1;
  ipf3[1] = 0;
  ipf3[2] = 0;
  ItFace it_face21_11_10 (2,1,periodic,n3,index,ic3,ipf3); 
  while (it_face21_11_10.next(if3)) {
    iface[IF3(if3)] = true;
    printf ("set %d %d %d\n",if3[0],if3[1],if3[2]); 
    ++count;
  }
  unit_assert(count == 1);
  unit_assert(iface[IF(1, 0, 0)]);


  count = 0;
  clear(iface,27);
  unit_func("ItFace(2,0,child(1,1),face(1,0))");
  ic3[0] = 1;
  ic3[1] = 1;
  ic3[2] = 0;
  ipf3[0] = 1;
  ipf3[1] = 0;
  ipf3[2] = 0;
  ItFace it_face20_11_10 (2,0,periodic,n3,index,ic3,ipf3); 
  while (it_face20_11_10.next(if3)) {
    iface[IF3(if3)] = true;
    printf ("set %d %d %d\n",if3[0],if3[1],if3[2]); 
    ++count;
  }
  unit_assert(count == 2);
  unit_assert(iface[IF(1, 0, 0)] && iface[IF(1,-1,0)]);

  // -------- 3D --------

  count = 0;
  clear(iface,27);
  unit_func("ItFace(3,0)");
  ItFace it_face30 (3,0,periodic,n3,index);
  while (it_face30.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 26);

  count = 0;
  clear(iface,27);
  ic3[0] = 1;
  ic3[1] = 0;
  ic3[2] = 1;
  unit_func("ItFace(3,0,child(1,0,1))");
  ItFace it_face30_101 (3,0,periodic,n3,index,ic3); 
  while (it_face30_101.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 7);
  unit_assert(iface[IF(0,-1,0)] && 
	      iface[IF(0,-1,1)] && 
	      iface[IF(0, 0,1)] &&
	      iface[IF(1,-1,0)] && 
	      iface[IF(1,-1,1)] && 
	      iface[IF(1, 0,0)] && 
	      iface[IF(1, 0,1)]);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(3,1)");
  ItFace it_face31 (3,1,periodic,n3,index); 
  while (it_face31.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 18);

  count = 0;
  clear(iface,27);
  ic3[0] = 0;
  ic3[1] = 0;
  ic3[2] = 1;
  unit_func("ItFace(3,1,child(0,0,1))");
  ItFace it_face31_001 (3,1,periodic,n3,index,ic3); 
  while (it_face31_001.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 6);
  unit_assert(iface[IF(-1,-1,0)] && 
	      iface[IF(-1, 0,1)] &&
	      iface[IF(-1, 0,0)] && 
	      iface[IF( 0,-1,0)] && 
	      iface[IF( 0,-1,1)] && 
	      iface[IF( 0, 0,1)]);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(3,2)");
  ItFace it_face32 (3,2,periodic,n3,index); 
  while (it_face32.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 6);

  count = 0;
  clear(iface,27);
  ic3[0] = 1;
  ic3[1] = 1;
  ic3[2] = 0;
  unit_func("ItFace(3,2,child_110)");
  ItFace it_face32_110 (3,2,periodic,n3,index,ic3); 
  while (it_face32_110.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 3);
  unit_assert(iface[IF( 0, 0,-1)] && 
	      iface[IF( 0, 1, 0)] &&
	      iface[IF( 1, 0, 0)]);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(3,1,child(0,1,1),face(0,0,1))");
  ic3[0] = 0;
  ic3[1] = 1;
  ic3[2] = 1;
  ipf3[0] = 0;
  ipf3[1] = 0;
  ipf3[2] = 1;
  ItFace it_face31_011_001 (3,1,periodic,n3,index,ic3,ipf3); 
  while (it_face31_011_001.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 3);
  unit_assert(iface[IF(0,0,1)] && iface[IF(0,-1,1)] &&
	      iface[IF(1,0,1)]);

  count = 0;
  clear(iface,27);
  unit_func("ItFace(3,0,child(1,0,0),face(1,1,0))");
  ic3[0] = 1;
  ic3[1] = 1;
  ic3[2] = 0;
  ipf3[0] = 1;
  ipf3[1] = -1;
  ipf3[2] = 0;
  ItFace it_face30_100_110 (3,0,periodic,n3,index,ic3,ipf3); 
  while (it_face30_100_110.next(if3)) {
    iface[IF3(if3)] = true;
    ++count;
  }
  unit_assert(count == 2);
  unit_assert(iface[IF(1,-1,0)] && iface[IF(1,-1,1)]);

  //--------------------------------------------------

  //--------------------------------------------------

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

