#ifndef NODE16_H
#define NODE16_H

#include <stdlib.h>

// Production rules for ordering_hilbert4[type] 
const int ordering_hilbert4[4][4] = 
  { {3,0,0,1},
    {1,1,2,0},
    {2,3,1,2},
    {0,2,3,3}};


class Tree16;

// WARNING: NUMBERING MUST BE SUCH THAT DIR AND (DIR+2) % 4 MUST BE
// OPPOSITE DIRECTIONS

enum face_type {
  R = 0,
  U = 1,
  L = 2,
  D = 3 };

enum corner_type {
  UL = 0,
  DL = 1,
  UR = 2,
  DR = 3 };

enum edge_type {
  LL = 0,
  CL = 1,
  CR = 2,
  RR = 3 };

face_type opposite(face_type face);

class Node16 {

 public:

  // Create a new leaf node
  Node16( int level_adjust = 0 );

  // Delete a node and all descedents
  ~Node16();

  // return the num'th child
  Node16 * child (int ix, int iy);

  // return the num'th neighbor
  Node16 * neighbor (face_type face);

  // make the two nodes neighbors.  friend function since either can be NULL
  friend void make_neighbors (Node16 * node_1, face_type face_1, Node16 * node_2);

  // get the child's cousin
  Node16 * cousin (face_type face, int ix, int iy);

  // return the parent
  Node16 * parent ();

  // Refine if any elements in the array are true and recurse
  // return the level
  int refine 
    (
     const int * level_array, 
     int nd0,  int nd1,
     int low0, int up0,  
     int low1, int up1,
     int level, 
     int max_level,
     bool is_full = true
     );

  // Perform a pass of trying to remove level-jumps 
  void normalize_pass(bool & refined_tree, bool is_full = true);

  // Perform a pass of trying to optimize uniformly-refined nodes
  void optimize_pass(bool & refined_tree, bool is_full = true);

  // Fill the image region with values
  void fill_image
    (
     float * image,
     int nd0,  int nd1,
     int low0, int up0,  
     int low1, int up1,
     int level,
     int num_levels,
     int line_width
     );

  // Return whether node has all children
  bool all_children () {
    return 
      ((child_[0][0]) && (child_[1][0]) && (child_[2][0]) && (child_[3][0]) &&
       (child_[0][1]) && (child_[1][1]) && (child_[2][1]) && (child_[3][1]) &&
       (child_[0][2]) && (child_[1][2]) && (child_[2][2]) && (child_[3][2]) &&
       (child_[0][3]) && (child_[1][3]) && (child_[2][3]) && (child_[3][3]));
  };

  // Return whether node has any children
  bool any_children () { 
    return 
      ((child_[0][0]) || (child_[1][0]) || (child_[2][0]) || (child_[3][0]) ||
       (child_[0][1]) || (child_[1][1]) || (child_[2][1]) || (child_[3][1]) ||
       (child_[0][2]) || (child_[1][2]) || (child_[2][2]) || (child_[3][2]) ||
       (child_[0][3]) || (child_[1][3]) || (child_[2][3]) || (child_[3][3]));
  };

  static int num_nodes() { return num_nodes_; };

 private:

  Node16 * child_[4][4];    // column-ordering
  Node16 * neighbor_[4]; // Right up left down
  Node16 * parent_;
  int level_adjust_;      // scale for optimizing uniformly refined nodes


  void create_children_();
  void update_children_();
  void delete_children_();

  void update_child_ (int ix, int iy);
  void create_child_ (int ix, int iy);

  static int num_nodes_;

};

#endif
