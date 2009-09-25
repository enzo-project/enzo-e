//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2009 James Bordner
 * Copyright (C) 2009 Laboratory for Computational Astrophysics
 * Copyright (C) 2009 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      test_memory.cpp
 * @brief     Program implementing unit tests for the Memory component
 * @author    James Bordner
 * @date      Thu Sep  3 16:37:24 PDT 2009
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include <stdio.h>
#include <stdlib.h>

#include <string>

#include "cello.h"
#include "performance_timer.hpp"
#include "error.hpp"
#include "memory.hpp"
#include "test.hpp"

main()
{

  unit_class ("Memory");
  unit_open();

  //----------------------------------------------------------------------
  // allocate()
  //----------------------------------------------------------------------

  unit_func("allocate");

  int i=0;
  double *f1;
  float  *f2;
  int    *f3;

  size_t size = 0;

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

  NEW_F1(size);
  NEW_F2(size);
  NEW_F3(size);

  unit_assert (Memory::bytes() == size);
  
  //----------------------------------------------------------------------
  // deallocate()
  //----------------------------------------------------------------------

  unit_func("deallocate");

  DEL_F1(size);
  DEL_F3(size);

  unit_assert (Memory::bytes() == size);

  DEL_F2(size);
  NEW_F1(size);
  NEW_F2(size);
  DEL_F1(size);
  DEL_F2(size);
  NEW_F3(size);
  NEW_F2(size);
  DEL_F3(size);
  DEL_F2(size);

  unit_assert (Memory::bytes() == size);

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
  

  unit_assert (Memory::bytes() == size);

  //----------------------------------------------------------------------
  // begin_group(), end_group()
  //----------------------------------------------------------------------

  size_t size_test_1 = 0;
  size_t size_test_2 = 0;

  // Group 1

  Memory::begin_group("Test 1");

  int handle_1 = Memory::current_handle();

  NEW_F1(size_test_1);
  NEW_F3(size_test_1);

  unit_assert (Memory::bytes(handle_1) == size_test_1);
  unit_assert (Memory::bytes() == size);
  
  DEL_F1(size_test_1);
  DEL_F3(size_test_1);

  unit_assert (Memory::bytes(handle_1) == size_test_1);
  unit_assert (Memory::bytes() == size);
  unit_assert (strcmp(Memory::current_group(),"Test 1") == 0);

  Memory::end_group("Test 1");

  // Group 1

  Memory::begin_group("Test 2");

  int handle_2 = Memory::current_handle();

  NEW_F2(size_test_2);
  NEW_F3(size_test_2);

  unit_assert (Memory::bytes(handle_1) == size_test_1);
  unit_assert (Memory::bytes(handle_2) == size_test_2);
  unit_assert (Memory::bytes() == size);
  unit_assert (strcmp(Memory::current_group(),"Test 2") == 0);
  
  Memory::end_group("Test 2");

  DEL_F2(size_test_2);
  DEL_F3(size_test_2);
 
  unit_assert (Memory::bytes(handle_1) == size_test_1);
  unit_assert (Memory::bytes(handle_2) == size_test_2);
  unit_assert (Memory::bytes() == size);
  unit_assert (strcmp(Memory::current_group(),"\0") == 0);

  // curr_group()
  // bytes()
  // available()
  // efficiency()
  // highest()
  // set_limit()
  // get_limit()
  // num_new()
  // num_delete()

  unit_close();
}
