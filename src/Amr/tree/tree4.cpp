
#include <stdio.h>
#include "cello.h"
#include "tree4.h"
#include "node4.h"

Tree4::Tree4()
  : root_(new Node4(0)),
    levels_(0)
{

}

void Tree4::refine(const bool * mask_array, int nd0, int nd1)
{
  levels_ = root_->refine(mask_array,nd0, nd1,0,nd0,0,nd1);
}

/// Print the tree
void Tree4::print()
{
  root_->print(0);
}

/// Create an hdf5 file of tree, assuming given source bitmap size
/// 
float * Tree4::create_image 
(int * nd0, int * nd1, int n0, int n1,int cell_width)
{
  *nd0 = (n0+1) + cell_width*n0;
  *nd1 = (n1+1) + cell_width*n1;
  float * image = new float [*nd0*(*nd1)];
  
  root_->fill_image(image,(*nd0),(*nd1),0,(*nd0)-1,0,(*nd1)-1,0,levels_);
  return image;
}

