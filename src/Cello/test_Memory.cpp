// See LICENSE_CELLO file for license and copyright information

/// @file      test_Memory.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Sep  3 16:37:24 PDT 2009
/// @brief     Program implementing unit tests for the Memory component
 
#include "main.hpp" 
#include "test.hpp"

#include "performance.hpp" /* for Timer */
#include "memory.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Memory");

#ifdef CONFIG_USE_MEMORY
  Memory * memory = Memory::instance();

  memory->reset();

  PARALLEL_PRINTF ("start\n"); fflush(stdout);

  //----------------------------------------------------------------------
  // allocate()
  //----------------------------------------------------------------------

  unit_func("allocate");

  int i=0;
  double *f1;
  float  *f2;
  int    *f3;
  int new_count = 0;
  int del_count = 0;

#define NEW(VAR,TYPE,SIZE,COUNT) \
  VAR = new TYPE[SIZE]; \
  COUNT += sizeof(TYPE[SIZE]); \
  for (i=0; i<SIZE; i++) VAR[i] = 17; \
  new_count++;

#define DEL(VAR,TYPE,SIZE,COUNT) \
  delete [] VAR; \
  COUNT -= sizeof(TYPE[SIZE]); \
  del_count++;

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
  const int size_alloc = 1000;
  for (int j=0; j<num_alloc; j++) {
    f1 = new double[size_alloc];
    new_count++;
    delete [] f1;
    del_count++;
  }
  timer.stop();
  PARALLEL_PRINTF ("new/delete  per sec = %g\n",num_alloc/timer.value());

  timer.clear();

  timer.start();
  for (int j=0; j<num_alloc; j++) {
    f1 = (double*)malloc(sizeof(double[size_alloc]));
    free(f1);
  }
  timer.stop();
  PARALLEL_PRINTF ("malloc/free per sec = %g\n",num_alloc/timer.value());

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

  const char * name_1;
  int handle_1 = memory->current_group(&name_1);

  unit_assert (strcmp(name_1,"Test_1") == 0);

  NEW_F1(size_test_1);
  NEW_F3(size_test_1);

  unit_assert (memory->bytes(handle_1) == size_test_1);
  printf ("%lld %d\n",memory->bytes()         , size_test_1 + size);
  unit_assert (memory->bytes()         == size_test_1 + size);
  
  DEL_F1(size_test_1);
  DEL_F3(size_test_1);

  unit_assert (memory->bytes(handle_1) == size_test_1);
  unit_assert (memory->bytes()         == size_test_1 + size);

  memory->end_group(group_test_1);

  // Group 1

  memory->begin_group(group_test_2);

  const char * name_2;
  int handle_2 = memory->current_group(&name_2);

  unit_assert (strcmp(name_2,"Test_2") == 0);

  NEW_F2(size_test_2);
  NEW_F3(size_test_2);

  unit_assert (memory->bytes(handle_1) == size_test_1);
  unit_assert (memory->bytes(handle_2) ==               size_test_2);
  unit_assert (memory->bytes()         == size_test_1 + size_test_2 + size);
  
  memory->end_group(group_test_2);

  const char * name_0;
  int handle_0 = memory->current_group(&name_0);

  unit_assert (strcmp(name_0,"\0") == 0);
  unit_assert (handle_0 == 0);

  DEL_F2(size_test_2);
  DEL_F3(size_test_2);

  unit_assert (memory->bytes(handle_1) == size_test_1);
  unit_assert (memory->bytes(handle_2) ==               size_test_2);
  unit_assert (memory->bytes()         == size_test_1 + size_test_2 + size);

  // available()

  memory->set_limit(1000000);
  memory->set_limit(10000,1);

  unit_func ("set_limit()");
  unit_assert(memory->limit()  == 1000000);
  unit_func ("limit()");
  unit_assert(memory->limit(1) == 10000);

  // efficiency()

  
  unit_func ("efficiency()");

  char * temp_0 = new char [10000];
  new_count++;
  unit_assert (fabs(memory->efficiency() - 0.01) < 1e-7);
  delete [] temp_0;
  del_count++;

  memory->begin_group(1);
  char * temp_1 = new char [1000];
  new_count++;
  unit_assert (fabs(memory->efficiency(1) - 0.1) < 1e-7);
  memory->end_group(1);

  delete [] temp_1;
  del_count++;

  // bytes_high()
  unit_func ("bytes_high()");

  unit_assert(memory->bytes_high() == 10000);
  unit_assert(memory->bytes_high(1) == 1000);

  // num_new()
  unit_func ("num_new()");
  unit_assert(memory->num_new() == new_count);

  // num_delete()
  unit_func ("num_delete()");
  unit_assert(memory->num_delete() == del_count);

  memory->print();
#else /* CONFIG_USE_MEMORY */
  unit_func("CONFIG_USE_MEMORY");
  unit_assert(true);
#endif/* CONFIG_USE_MEMORY */

  unit_finalize();

  PARALLEL_EXIT;

}

PARALLEL_MAIN_END

