// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      test_memory.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Sep  3 16:37:24 PDT 2009
/// @brief     Program implementing unit tests for the Memory component
 
#include <stdio.h>
#include <stdlib.h>

#include <string>

#include "cello.h"

#include "error.hpp"
#include "performance.hpp"
#include "memory.hpp"
#include "test.hpp"

int main(int argc, char ** argv)
{

  unit_init();

  Memory * memory = Memory::instance();

  printf ("start\n"); fflush(stdout);

  unit_class ("Memory");

  //----------------------------------------------------------------------
  // allocate()
  //----------------------------------------------------------------------

  unit_func("allocate");

  int i=0;
  double *f1;
  float  *f2;
  int    *f3;

#define NEW(VAR,TYPE,SIZE,COUNT)			\
  VAR = new TYPE[SIZE]; \
  COUNT += sizeof(TYPE[SIZE]);			\
  for (i=0; i<SIZE; i++) VAR[i] = 17;

#define DEL(VAR,TYPE,SIZE,COUNT)			\
  delete [] VAR; \
  COUNT -= sizeof(TYPE[SIZE]);

#define NEW_F1(SIZE) NEW(f1,double,10,SIZE);
#define NEW_F2(SIZE) NEW(f2,float, 17,SIZE);
#define NEW_F3(SIZE) NEW(f3,int,   25,SIZE);
#define DEL_F1(SIZE) DEL(f1,double,10,SIZE);
#define DEL_F2(SIZE) DEL(f2,float, 17,SIZE);
#define DEL_F3(SIZE) DEL(f3,int,   25,SIZE);

  unsigned size = 0;

  NEW_F1(size);
  NEW_F2(size);
  NEW_F3(size);

  unit_assert (memory->bytes() == size);
  
  //----------------------------------------------------------------------
  // deallocate()
  //----------------------------------------------------------------------

  unit_func("deallocate");

  DEL_F1(size);
  DEL_F3(size);

  unit_assert (memory->bytes() == size);

  DEL_F2(size);
  NEW_F1(size);
  NEW_F2(size);
  DEL_F1(size);
  DEL_F2(size);
  NEW_F3(size);
  NEW_F2(size);
  DEL_F3(size);
  DEL_F2(size);

  unit_assert (memory->bytes() == size);

  Timer timer;
  timer.start();
  const int num_alloc = 10000;
  const int size_alloc = 1000000;
  for (int j=0; j<num_alloc; j++) {
    f1 = new double[size_alloc];
    delete [] f1;
  }
  timer.stop();
  printf ("alloc/dealloc per sec = %g\n",num_alloc/timer.value());

  timer.clear();
  timer.start();
  for (int j=0; j<num_alloc; j++) {
    f1 = (double*)malloc(sizeof(double[size_alloc]));
    free(f1);
  }
  timer.stop();
  printf ("  new/delete  per sec = %g\n",num_alloc/timer.value());
  

  unit_assert (memory->bytes() == size);

  //----------------------------------------------------------------------
  // begin_group(), end_group()
  //----------------------------------------------------------------------

  unsigned size_test_1 = 0;
  unsigned size_test_2 = 0;

  // Group 1

  unsigned group_test_1 = 1;
  unsigned group_test_2 = 2;

  memory->new_group (group_test_1,"Test_1");
  memory->new_group (group_test_2,"Test_2");

  memory->begin_group(group_test_1);

  unit_assert (strcmp(memory->current_group(),"Test_1") == 0);

  int handle_1 = memory->current_handle();

  NEW_F1(size_test_1);
  NEW_F3(size_test_1);

  unit_assert (memory->bytes(handle_1) == size_test_1);
  unit_assert (memory->bytes() == size + size_test_1);
  
  DEL_F1(size_test_1);
  DEL_F3(size_test_1);

  unit_assert (memory->bytes(handle_1) == size_test_1);
  unit_assert (memory->bytes() == size);

  memory->end_group(group_test_1);

  // Group 1

  memory->begin_group(group_test_2);

  unit_assert (strcmp(memory->current_group(),"Test_2") == 0);

  int handle_2 = memory->current_handle();

  NEW_F2(size_test_2);
  NEW_F3(size_test_2);

  unit_assert (memory->bytes(handle_1) == size_test_1);
  unit_assert (memory->bytes(handle_2) == size_test_2);
  unit_assert (memory->bytes() == size + size_test_1 + size_test_2);
  
  memory->end_group(group_test_2);

  unit_assert (strcmp(memory->current_group(),"\0") == 0);

  DEL_F2(size_test_2);
  DEL_F3(size_test_2);
 
  unit_assert (memory->bytes(handle_1) == size_test_1);
  unit_assert (memory->bytes(handle_2) == size_test_2);
  unit_assert (memory->bytes() == size);


  // curr_group()
  unit_func ("curr_group()");
  unit_assert(0);  //FAILS

  // bytes()
  unit_func ("bytes()");
  unit_assert(0); //FAILS

  // available()
  unit_func ("available()");
  unit_assert(0); //FAILS

  // efficiency()
  unit_func ("efficiency()");
  unit_assert(0); //FAILS

  // highest()
  unit_func ("highest()");
  unit_assert(0); //FAILS

  // set_limit()
  unit_func ("set_limit()");
  unit_assert(0); //FAILS

  // get_limit()
  unit_func ("get_limit()");
  unit_assert(0); //FAILS

  // num_new()
  unit_func ("num_new()");
  unit_assert(0); //FAILS

  // num_delete()
  unit_func ("num_delete()");
  unit_assert(0); //FAILS

  memory->print();

  unit_finalize();
}
