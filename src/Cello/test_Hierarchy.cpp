// See LICENSE_CELLO file for license and copyright information

/// @file     test_Hierarchy.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Hierarchy class

#include "main.hpp"
#include "test.hpp"
#include "mesh.hpp"

#include "mesh_charm.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;
  unit_init(0,1);

  unit_class("Hierarchy");

  unit_func("Hierarchy");
  Factory * factory = new Factory;
  const int dimensions = 3;
  const int refinement = 2;
  Hierarchy * hierarchy = factory->create_hierarchy(dimensions, refinement);
  unit_assert(hierarchy != NULL);

  FieldDescr field_descr;

  hierarchy->create_root_patch(&field_descr,12,12,12,3,3,3);
  unit_assert(hierarchy->patch(0)!=NULL);

  // Extents

  hierarchy->set_lower(-1.0, -2.0, -3.0);
  double xm,ym,zm;
  hierarchy->lower(&xm,&ym,&zm);
  unit_func("lower");
  unit_assert(xm == -1.0 && ym == -2.0 && zm == -3.0);

  hierarchy->set_upper( 1.0,  2.0,  3.0);
  double xp,yp,zp;
  hierarchy->upper(&xp,&yp,&zp);
  unit_func("upper");
  unit_assert(xp ==  1.0 && yp ==  2.0 && zp ==  3.0);

  
  // // Write header data
  // hierarchy->open(file,"hierarchy_header.out","w");
  // hierarchy->write(file,file_content_header);
  // hierarchy->close(file);

  // // Read header data
  // hierarchy_read->open(file,"hierarchy_header.out","r");
  // hierarchy_read->read(file,file_content_header);
  // hierarchy_read->close(file);
    
  unit_finalize();

  delete hierarchy;
  delete factory;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END
