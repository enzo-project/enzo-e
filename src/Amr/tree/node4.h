#include <stdlib.h>

// Production rules for ordering_hilbert4[type] 
const int ordering_hilbert4[4][4] = 
  { {3,0,0,1},
    {1,1,2,0},
    {2,3,1,2},
    {0,2,3,3}};


class Tree4;

enum face_type {
  face_R = 0,
  face_U = 1,
  face_L = 2,
  face_D = 3 };

enum corner_type {
  corner_UL = 0,
  corner_DL = 1,
  corner_UR = 2,
  corner_DR = 3 };

class Node4 {

 public:

  // Create a new leaf node
  Node4();

  // return the num'th child
  Node4 * child (int num);

  // return the num'th neighbor
  Node4 * neighbor (int num);

  // return the parent
  Node4 * parent ();

  // Refine if any elements in the array are true and recurse
  // return the level
  int refine 
    (
     const bool * array, 
     int nd0,  int nd1,
     int low0, int up0,  
     int low1, int up1
     );

  // Perform a pass of trying to remove level-jumps 
  bool normalize_pass();

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

  bool is_leaf() { return child_[0] == NULL; };


 private:

  Node4 * child_[4];    // column-ordering
  Node4 * neighbor_[4]; // Right up left down
  Node4 * parent_;

  void create_children_();

};

