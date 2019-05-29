// See LICENSE_CELLO file for license and copyright information

/// @file     test_Node.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-09
/// @brief    Test program for the Node class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Sync");

  Sync sync_1;
  Sync sync_2(10);

  unit_func("value");

  unit_assert(sync_1.value() == 0);
  unit_assert(sync_2.value() == 0);

  unit_func("stop");
  unit_assert(sync_1.stop() == 0);
  unit_assert(sync_2.stop() == 10);

  unit_func("set_stop");
  sync_1.set_stop(5);
  unit_assert(sync_1.stop() == 5);
  
  unit_func("next");
  unit_assert(sync_1.next() == false);
  unit_assert(sync_1.value() == 1);
  unit_assert(sync_1.is_done() == false);
  unit_assert(sync_1.next() == false);
  unit_assert(sync_1.value() == 2);
  unit_assert(sync_1.is_done() == false);
  unit_assert(sync_1.next() == false);
  unit_assert(sync_1.value() == 3);
  unit_assert(sync_1.is_done() == false);
  unit_assert(sync_1.next() == false);
  unit_assert(sync_1.value() == 4);
  unit_assert(sync_1.is_done() == false);
  unit_assert(sync_1.next() == true);
  unit_assert(sync_1.value() == 5);
  unit_assert(sync_1.is_done() == true);

  sync_1.clear();
  unit_assert(sync_1.is_done() == false);
  unit_assert(sync_1.value() == 0);
  unit_assert(sync_1.stop() == 5);
  sync_1.reset();
  unit_assert(sync_1.value() == 0);
  unit_assert(sync_1.stop() == 0);
  
  //   Sync (int index_stop = 0) throw()
  //  inline bool next () throw()
  //  inline bool is_done () const throw()
  //  inline void set_stop (int stop) throw ()
  //  inline void inc_stop (int increment) throw ()
  //  inline int value () const
  //  inline int stop () const throw ()
  //  inline void reset () throw () 
  //  inline int operator -- () 
  //  inline int operator ++ () 
  //  inline int operator -= (int count) 
  //  inline int operator += (int count) 

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

