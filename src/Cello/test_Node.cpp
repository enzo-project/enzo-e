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

  unit_class("Node");

  Node * root = new Node;

  //--------------------------------------------------

  unit_func("Node");

  unit_assert (root != NULL);

  unit_func("sizeof(Node) <= 48");
  unit_assert (sizeof(*root) <= 48);
  printf ("%lu\n",sizeof(*root));

  unit_func("is_leaf()");
  unit_assert (root->is_leaf() == true);

  //--------------------------------------------------

  unit_func("set_data()");

  double data = 10.0;
  root->set_data(&data);

  unit_assert (true);
  
  //--------------------------------------------------

  unit_func("data()");

  unit_assert (*((double *)(root->data())) == 10.0);

  //--------------------------------------------------

  unit_func("refine()");

  int root_children = 4;

  root->refine(root_children);

  for (int i=0; i<root_children; i++) 
    unit_assert (root->child(i) != NULL);

  Node * child = root->child(root_children-1);

  unit_func("is_leaf()");
  unit_assert (root->is_leaf()  == false);
  unit_assert (child->is_leaf() == true);

  unit_func("refine()");

  int child_children = 64;

  child->refine(child_children);

  for (int i=0; i<child_children; i++) 
    unit_assert (child->child(i) != NULL);

  //--------------------------------------------------

  unit_func ("is_leaf()");

  unit_assert (root->is_leaf() == false);
  unit_assert (child->is_leaf() == false);
  unit_assert (child->child(child_children-1)->is_leaf() == true);

  //--------------------------------------------------

  unit_func("coarsen()");

  child->coarsen(child_children);

  unit_assert(child->child(0) == NULL);
  unit_assert(root->child(0)  != NULL);

  root->coarsen(root_children);

  unit_assert(root->child(0)  == NULL);

  //--------------------------------------------------

  unit_func("~Node()");
  delete (root);
  unit_assert (true);

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

