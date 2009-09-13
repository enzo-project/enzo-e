#include <stdlib.h>

// Production rules for ordering_hilbert4[type] 
const int ordering_hilbert4[4][4] = 
  { {3,0,0,1},
    {1,1,2,0},
    {2,3,1,2},
    {0,2,3,3}};

class Tree4;

class Node4 {

 public:

  // Create a new leaf node
  Node4(int type);

  // return the num'th child
  Node4 * child (int num);

  // Print the tree to stdout
  void print(int level);
  
  // Refine if any elements in the array are true and recurse
  // return the level
  int refine 
    (
     const bool * array, 
     int nd0,  int nd1,
     int low0, int up0,  
     int low1, int up1
     );

  // Fill the image region with values
  void fill_image
    (
     float * image,
     int nd0,  int nd1,
     int low0, int up0,  
     int low1, int up1,
     int level,
     int num_levels
     );

 private:

  // types:
  // 0 cup open up
  // 1 cup open right
  // 2 cup open down
  // 3 cup open left

  int type_; // 0,1,2,3
  Node4 * children_[4]; // column-ordering

};

