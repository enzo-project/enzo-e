// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Parse.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-06
/// @brief    Test program for reading in parameters then displaying them

#include "parameters.hpp"

int main(int argc, char **argv)
{
  cello_parameters_read(stdin);
  cello_parameters_print();
  return 0;
}
