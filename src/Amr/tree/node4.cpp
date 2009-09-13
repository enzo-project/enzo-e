#include <stdio.h>
#include "cello.h"
#include "node4.h"

const bool debug = false;

Node4::Node4(int type) 
  : type_(type) 
{ 
  children_[0] = NULL;
  children_[1] = NULL;
  children_[2] = NULL;
  children_[3] = NULL;
}

Node4 * Node4::child (int num) 
{ 
  return children_[num]; 
}

// Create 4 empty child nodes

int Node4::refine 
(
 const bool * mask_array, 
 int nd0,  int nd1,
 int low0, int up0,  
 int low1, int up1
)
{
  int level = 0;
  if (up0-1 > low0 && up1-1 > low1) {

    if (debug) printf ("refine (%d %d  %d:%d %d:%d)\n",nd0,nd1,low0,up0,low1,up1);

    // determine whether to refine the node

    bool refine_node = false;

    for (int i1=low1; i1<up1; i1++) {
      for (int i0=low0; i0<up0; i0++) {
	int i = i0 + nd0 * i1;
	if (mask_array[i]) {
	  refine_node = true;
	  break;
	}
      }
    }

    // refine the node if needed

    if (refine_node) {
      int mid0 = (up0 + low0)/2;
      int mid1 = (up1 + low1)/2;

      
      // top-left
      children_[0] = new Node4(ordering_hilbert4[type_][0]);
      int l0 = children_[0]->refine (mask_array,nd0,nd1,low0,mid0,low1,mid1);

      // bottom-left
      children_[1] = new Node4(ordering_hilbert4[type_][1]);
      int l1 = children_[1]->refine (mask_array,nd0,nd1,mid0,up0,low1,mid1);

      // bottom-right
      children_[2] = new Node4(ordering_hilbert4[type_][2]);
      int l2 = children_[2]->refine (mask_array,nd0,nd1,low0,mid0,mid1,up1);

      // top-right
      children_[3] = new Node4(ordering_hilbert4[type_][3]);
      int l3 = children_[3]->refine (mask_array,nd0,nd1,mid0,up0,mid1,up1);

      level = (l0 > level) ? l0 : level;
      level = (l1 > level) ? l1 : level;
      level = (l2 > level) ? l2 : level;
      level = (l3 > level) ? l3 : level;

      ++ level;

    } // if refine node

  } // if not at bottom of recursion
  return level;
}

void Node4::print(int level)
{
  printf ("Level %d  type %d\n",level,type_);
  if (children_[0] != NULL) children_[0]->print (level+1);
  if (children_[1] != NULL) children_[1]->print (level+1);
  if (children_[2] != NULL) children_[2]->print (level+1);
  if (children_[3] != NULL) children_[3]->print (level+1);
  
}

  // Fill the image region with values
void Node4::fill_image
(
 float * image,
 int nd0,  int nd1,
 int low0, int high0,  
 int low1, int high1,
 int level,
 int num_levels
 )
{
  int i0,i1,i;
  // Fill interior
  for (i1=low1; i1<=high1; i1++) {
    for (i0=low0; i0<=high0; i0++) {
      i = i0 + nd0 * i1;
      image[i] = num_levels - level + 1;
    }
  }

  // Draw border
  for (i0=low0; i0<=high0; i0++) {
    i1 = low1;
    image[i0 + nd0*i1] = 0;
    i1 = high1;
    image[i0 + nd0*i1] = 0;
  }

  for (i1=low1; i1<=high1; i1++) {
    i0 = low0;
    image[i0 + nd0*i1] = 0;
    i0 = high0;
    image[i0 + nd0*i1] = 0;
  }
    

  // Recurse

  int mid0 = (high0 + low0)/2;
  int mid1 = (high1 + low1)/2;

  if (children_[0] != NULL) {
    children_[0]->fill_image (image,nd0,nd1,low0,mid0, low1,mid1, level + 1, num_levels);
  }
  if (children_[1] != NULL) {
    children_[1]->fill_image (image,nd0,nd1,mid0,high0,low1,mid1, level + 1, num_levels);
  }
  if (children_[2] != NULL) {
    children_[2]->fill_image (image,nd0,nd1,low0,mid0, mid1,high1,level + 1, num_levels);
  }
  if (children_[3] != NULL) {
    children_[3]->fill_image (image,nd0,nd1,mid0,high0,mid1,high1,level + 1, num_levels);
  }
}

