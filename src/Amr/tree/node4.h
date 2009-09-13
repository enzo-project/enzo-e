#include <stdlib.h>

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

  /// Create a new leaf node
  Node4(int type) : type_(type) 
    { 
      children_[0] = NULL;
      children_[1] = NULL;
      children_[2] = NULL;
      children_[3] = NULL;
    };

    /// return the num'th child
    Node4 * child (int num) { return children_[num]; }

    /// Create 4 empty child nodes
    void refine() {
      children_[0] = new Node4(ordering_hilbert4[type_][0]);
      children_[1] = new Node4(ordering_hilbert4[type_][1]);
      children_[2] = new Node4(ordering_hilbert4[type_][2]);
      children_[3] = new Node4(ordering_hilbert4[type_][3]);
    }

 private:

  // types:
  // 0 cup open up
  // 1 cup open right
  // 2 cup open down
  // 3 cup open left

  int type_; // 0,1,2,3
  Node4 * children_[4]; // column-ordering

};

