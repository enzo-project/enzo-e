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
#include <string>

#include "error.hpp"
#include "memory.hpp"
#include "test.hpp"

main()
{

  unit_class ("Memory");
  unit_open();

  // allocate()

  unit_func("allocate");

  int i=0;
  double *f1;
  float  *f2;
  int    *f3;

  size_t size = 0;

  f1 = new double[10]; 
  size += sizeof(double[10]);
  for (i=0; i<10; i++) f1[i] = 7.0;

  f2 = new float[17];  
  size += sizeof(float[17]);
  for (i=0; i<17; i++) f2[i] = 8.0;

  f3 = new int[25];    
  size += sizeof(int[25]);
  for (i=0; i<25; i++) f3[i] = 9;

  unit_assert (Memory::current() == size);
  
  // deallocate()

  unit_func("deallocate");

  delete [] f1; 
  size -= sizeof(double[10]);
  delete [] f3; 
  size -= sizeof(int[25]);
  unit_assert (Memory::current() == size);

  delete [] f2; 
  size -= sizeof(float[17]);

  f1 = new double[10];
  size += sizeof(double[10]);
  for (i=0; i<10; i++) f1[i] = 7.0;

  f2 = new float; 
  size += sizeof(float);
  *f2 = 10;

  delete f1; 
  size -= sizeof(double[10]);

  delete f2; 
  size -= sizeof(float);

  f3 = new int[25]; 
  size += sizeof(int[25]);
  for (i=0; i<25; i++) f3[i] = 9;

  f2 = new float[17];
  size += sizeof(float[17]);
  for (i=0; i<17; i++) f2[i] = 8;
  
  delete [] f3; size -=  sizeof(int[25]);

  unit_assert (Memory::current() == size);

  // begin_group(), end_group()
  // curr_group()
  // current()
  // available()
  // efficiency()
  // highest()
  // set_limit()
  // get_limit()
  // num_new()
  // num_delete()

  unit_close();
}
