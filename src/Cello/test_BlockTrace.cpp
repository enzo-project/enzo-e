// See LICENSE_CELLO file for license and copyright information

/// @file     test_BlockTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-04-21
/// @brief    Test program for the BlockTrace class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("BlockTrace");

  // octree rooted at array element (3,2,1)
  Index index (2,2,4);

  int index_min[3] = {2,2,4};
  int index_max[3] = {4,5,6};
  BlockTrace * block_trace = new BlockTrace (3,index_min,index_max);

  // BlockTrace object was created
  unit_func ("BlockTrace()");
  unit_assert (block_trace != nullptr);

  // Home index is that associated with index_min[]
  unit_func("home()");
  unit_assert (block_trace->home() == index);
  
  // BlockTrace object root Index is correct
  unit_func ("top()");
  unit_assert (block_trace->top() == index);

  // traverse to child 0 
  unit_func ("next()");
  unit_assert ( ! block_trace->next(false) );

  // check with expected
  int level = 1;
  index.set_level(level);
  index.set_child(level,0,0,0);
  unit_assert (block_trace->top() == index);

  // traverse to child 00 
  unit_assert ( ! block_trace->next(false) );
  // check
  ++level;
  index.set_level(level);
  index.set_child(level,0,0,0);
  unit_assert (block_trace->top() == index);

  // check with expected
  // traverse to child 01 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,1,0,0);
  unit_assert (block_trace->top() == index);
  
  // traverse to child 02 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,0,1,0);
  unit_assert (block_trace->top() == index);

  // traverse to child 03 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,1,1,0);
  unit_assert (block_trace->top() == index);

  // traverse to child 04 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,0,0,1);
  unit_assert (block_trace->top() == index);

  // traverse to child 05 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,1,0,1);
  unit_assert (block_trace->top() == index);

  // traverse to child 06 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,0,1,1);
  unit_assert (block_trace->top() == index);

  // traverse to child 07 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,1,1,1);
  unit_assert (block_trace->top() == index);

  // traverse to child 1 
  unit_assert ( ! block_trace->next(true) );
  --level;
  index.set_level(level);
  index.set_child(level,1,0,0);
  unit_assert (block_trace->top() == index);

  // traverse to child 2 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,0,1,0);
  unit_assert (block_trace->top() == index);

  // traverse to child 3 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,1,1,0);
  unit_assert (block_trace->top() == index);

  // traverse to child 4 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,0,0,1);
  unit_assert (block_trace->top() == index);

  // traverse to child 5 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,1,0,1);
  unit_assert (block_trace->top() == index);

  // traverse to child 6 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,0,1,1);
  unit_assert (block_trace->top() == index);
  
  // traverse to child 7 
  unit_assert ( ! block_trace->next(true) );
  index.set_child(level,1,1,1);
  unit_assert (block_trace->top() == index);

  //--------------------

  unit_assert ( ! block_trace->next(true) );
  index.set_array(3,2,4);
  index.set_level(0);
  block_trace->print("last");
  unit_assert (block_trace->top() == index);
  
  // // traverse to root
  // unit_assert ( block_trace->next(true) );
  
  //--------------------------------------------------

  delete block_trace;

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

