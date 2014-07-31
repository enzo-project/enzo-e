// See LICENSE_CELLO file for license and copyright information

/// @file     test_ItField.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the ItField class

#include "main.hpp"
#include "test.hpp"

#include "field.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  //--------------------------------------------------
  unit_class("ItFieldRange");
  //--------------------------------------------------

  unit_func ("it_field_range");

  ItField * it_field_count = new ItFieldRange(10);

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

  ItField * it_field_range = new ItFieldRange (3,9);

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
  unit_class("ItFieldList");
  //--------------------------------------------------

  ItField * it_field_list = new ItFieldList();

  unit_func ("it_field_list");

  const int values[6] = {15,12,-17,2,2,0};

  ItFieldList * itcast = dynamic_cast<ItFieldList *>(it_field_list);
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
  unit_class("ItFieldRange");
  //--------------------------------------------------

  unit_finalize();
  
  exit_();
}

PARALLEL_MAIN_END

