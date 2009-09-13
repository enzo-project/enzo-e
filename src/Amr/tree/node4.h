
// Production rules for ordering_hilbert4[type] 
const int ordering_hilbert4[4][4] = 
  { {3,0,0,1},
    {1,1,2,0},
    {2,3,1,2},
    {0,2,3,3}};

class Tree4;

class Node4 {

  friend class Tree4;

 protected:

  Node4();

  Node4 * children[4]; // column-ordering

 private:

  // types:
  // 0 cup open up
  // 1 cup open right
  // 2 cup open down
  // 3 cup open left

  int type_; // 0,1,2,3

};

