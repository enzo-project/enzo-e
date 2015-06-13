// See LICENSE_CELLO file for license and copyright information

/// @file     test_Refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Refresh class

#include "main.hpp"
#include "test.hpp"

#include "problem.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Refresh");

  Refresh * refresh = new Refresh (3,2);

  unit_assert (refresh != NULL);

  //--------------------------------------------------

  unit_func ("field_ghost_depth()");
  unit_assert(refresh->ghost_depth() == 3);
  unit_func ("min_face_rank()");
  unit_assert(refresh->min_face_rank() == 2);

  unit_func ("add_field()");
  refresh->add_field (12);
  refresh->add_field (9);
  refresh->add_field (-2);
  unit_assert (refresh->field(12));
  unit_assert (! refresh->field(-8));
  unit_assert (! refresh->field(4));
  unit_assert (refresh->field(9));
  unit_assert (! refresh->field(0));
  unit_assert (refresh->field(-2));

  //--------------------------------------------------

  delete refresh;

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

