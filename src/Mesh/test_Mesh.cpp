// $Id: test_mesh.cpp 1302 2010-03-17 00:16:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Mesh class
 
#include <stdio.h>
#include <string>

#include "cello.hpp"
#include "test.hpp"
#include "mesh.hpp"

int main(int argc, char ** argv)
{
  unit_init();
  unit_class ("Mesh");
  unit_func("Mesh");
  Mesh mesh;
  unit_assert(false);
  unit_finalize();
}
