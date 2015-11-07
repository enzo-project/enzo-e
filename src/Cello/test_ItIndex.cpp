// See LICENSE_CELLO file for license and copyright information

/// @file     test_ItIndex.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the ItIndex class

#include "main.hpp"
#include "test.hpp"

#include "data.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  //--------------------------------------------------
  unit_class("ItIndexRange");
  //--------------------------------------------------

  unit_func ("it_field_range");

  ItIndex * it_field_count = new ItIndexRange(10);

  unit_assert (it_field_count != NULL);

  int s=0;
  for (it_field_count->first(); ! it_field_count->done(); it_field_count->next()) {
    s += it_field_count->value();
  }
  unit_assert (s == (9 * 10) / 2);

  s = 0;
  for (it_field_count->first(); ! it_field_count->done(); it_field_count->next()) {
    s += it_field_count->value();
  }
  unit_assert (s == (9 * 10) / 2);


  unit_func ("it_field_count");

  ItIndex * it_field_range = new ItIndexRange (3,9);

  s = 0;
  for (it_field_range->first(); ! it_field_range->done(); it_field_range->next()) {
    s += it_field_range->value();
  }

  unit_assert (s == 3+4+5+6+7+8+9);
  
  s = 0;
  for (it_field_range->first(); ! it_field_range->done(); it_field_range->next()) {
    s += it_field_range->value();
  }

  unit_assert (s == 3+4+5+6+7+8+9);
  
  //--------------------------------------------------
  unit_class("ItIndexList");
  //--------------------------------------------------

  ItIndex * it_field_list = new ItIndexList();

  unit_func ("it_field_list");

  const int values[6] = {15,12,-17,2,2,0};

  ItIndexList * itcast = dynamic_cast<ItIndexList *>(it_field_list);
  if ( itcast != NULL) {
    itcast->append(values[0]);
    itcast->append(values[1]);
    itcast->append(values[2]);
    itcast->append(values[3]);
    itcast->append(values[4]);
  }

  int i=0;
  for (it_field_list->first();
       ! it_field_list->done();
       it_field_list->next()) {
    unit_assert(it_field_list->value() == values[i]);
    ++i;
  }

  i=0;
  for (it_field_list->first();
       ! it_field_list->done();
       it_field_list->next()) {
    unit_assert(it_field_list->value() == values[i]);
    ++i;
  }
  
  //--------------------------------------------------
  unit_class("ItIndexRange");
  //--------------------------------------------------

  unit_finalize();
  
  exit_();
}

PARALLEL_MAIN_END

