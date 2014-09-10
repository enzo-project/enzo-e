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
  memory->set_active(true); \
  VAR = new TYPE[SIZE]; \
  COUNT += sizeof(TYPE[SIZE]); \
  for (i=0; i<SIZE; i++) VAR[i] = 17; \
  new_count++; \
  memory->set_active(false);

#define DEL(VAR,TYPE,SIZE,COUNT) \
  memory->set_active(true); \
  delete [] VAR; \
  COUNT -= sizeof(TYPE[SIZE]); \
  del_count++; \
  memory->set_active(false);

#define NEW_F1(SIZE) NEW(f1,double,10,SIZE);
#define NEW_F2(SIZE) NEW(f2,float, 17,SIZE);
#define NEW_F3(SIZE) NEW(f3,int,   25,SIZE);
#define DEL_F1(SIZE) DEL(f1,double,10,SIZE);
#define DEL_F2(SIZE) DEL(f2,float, 17,SIZE);
#define DEL_F3(SIZE) DEL(f3,int,   25,SIZE);

  unsigned bytes = 0;

  NEW_F1(bytes);
  NEW_F2(bytes);
  NEW_F3(bytes);

  unit_assert (memory->bytes() == bytes);
  
  //----------------------------------------------------------------------
  // deallocate()
  //----------------------------------------------------------------------

  unit_func("deallocate");

  DEL_F1(bytes);
  DEL_F3(bytes);

  unit_assert (memory->bytes() == bytes);

  DEL_F2(bytes);
  NEW_F1(bytes);
  NEW_F2(bytes);
  DEL_F1(bytes);
  DEL_F2(bytes);
  NEW_F3(bytes);
  NEW_F2(bytes);
  DEL_F3(bytes);
  DEL_F2(bytes);

  unit_assert (memory->bytes() == bytes);

  Timer timer;
  timer.start();
  const int num_alloc = 10000;
  const int size_alloc = 1000;
  for (int j=0; j<num_alloc; j++) {
    memory->set_active(true);
    f1 = new double[size_alloc];
    memory->set_active(false);
    new_count++;
    memory->set_active(true);
    delete [] f1;
    memory->set_active(false);
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

  unit_assert (memory->bytes() == bytes);

  //----------------------------------------------------------------------
  // set_group()
  //----------------------------------------------------------------------

  unsigned size_test_1 = 0;
  unsigned size_test_2 = 0;

  // Group 1

  memory->new_group ("Test_1");
  memory->new_group ("Test_2");

  memory->set_group("Test_1");

  std::string name_1 = memory->group();

  unit_assert (name_1 == "Test_1");

  NEW_F1(size_test_1);
  NEW_F3(size_test_1);

  long long bytes_1;
  bytes_1 = memory->bytes("Test_1");
  unit_assert (bytes_1 == size_test_1);
  unit_assert (memory->bytes()         == size_test_1 + bytes);
  
  DEL_F1(size_test_1);
  DEL_F3(size_test_1);

  unit_assert (memory->bytes("Test_1") == size_test_1);
  unit_assert (memory->bytes()         == size_test_1 + bytes);

  memory->set_group("Test_2");

  std::string name_2 = memory->group();

  unit_assert (name_2 == "Test_2");

  NEW_F2(size_test_2);
  NEW_F3(size_test_2);

  unit_assert (memory->bytes("Test_1") == size_test_1);
  unit_assert (memory->bytes("Test_2") ==               size_test_2);
  unit_assert (memory->bytes()         == size_test_1 + size_test_2 + bytes);
  
  DEL_F2(size_test_2);
  DEL_F3(size_test_2);

  unit_assert (memory->bytes("Test_1") == size_test_1);
  unit_assert (memory->bytes("Test_2") ==               size_test_2);
  unit_assert (memory->bytes()         == size_test_1 + size_test_2 + bytes);

  // available()

  memory->set_bytes_limit(1000000);
  memory->set_bytes_limit(10000,"Test_1");

  unit_func ("set_bytes_limit()");
  unit_assert(memory->bytes_limit()  == 1000000);
  unit_func ("bytes_limit()");
  unit_assert(memory->bytes_limit("Test_1") == 10000);

  // efficiency()

  
  unit_func ("efficiency()");

  memory->set_active(true);
  char * temp_0 = new char [10000];
  memory->set_active(false);

  new_count++;
  unit_assert (fabs(memory->efficiency() - 0.01) < 1e-7);
  memory->set_active(true);
  delete [] temp_0;
  memory->set_active(false);
  del_count++;

  memory->set_group("Test_1");
  memory->set_active(true);
  char * temp_1 = new char [1000];
  memory->set_active(false);
  new_count++;
  unit_assert (fabs(memory->efficiency("Test_1") - 0.1) < 1e-7);

  memory->set_active(true);
  delete [] temp_1;
  memory->set_active(false);
  del_count++;

  // bytes_high()
  unit_func ("bytes_high()");

  unit_assert(memory->bytes_high() == 10000);
  unit_assert(memory->bytes_high("Test_1") == 1000);

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

  exit_();

}

PARALLEL_MAIN_END

