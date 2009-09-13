
#include "cello.h"
#include "tree4.h"
#include "node4.h"

// create a tree refined to the array non-zeros
// assume array width and height are powers of 2
Tree4::Tree4(const bool *array, int width, int height)
  : root_(new Node4(0)),
    array_(array),
    width_(width),
    height_(height)
{

}

/// Print the tree
void Tree4::print()
{
}
