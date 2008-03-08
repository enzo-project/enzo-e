/** 
 *********************************************************************
 *
 * @file      test_hydro.cpp
 * @brief     Program implementing unit tests for hydrodynamics
 * @author    James Bordner
 * @date      Fri Mar  7 17:11:14 PST 2008
 *
 *********************************************************************
 */
 
 
#include <stdio.h>

#include "scalar.hpp"
#include "array.hpp"
#include "unit.hpp"

main()
{
  unit_class ("Array");
  unit_open();
  unit_class_size(Array);

}
