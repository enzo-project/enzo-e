#define MAX_LEVELS 80

class Node4;

class Tree4 {

 public:

  // create a tree refined to the array non-zeros
  // assume array width and height are powers of 2
  Tree4();

  // Refine down to array
  void refine
    (const bool * array, 
     int nd0, int nd1, 
     int max_level, 
     bool is_full = true
     );

  // print levels
  void print_levels();

  // Refine nodes to remove level jumps
  void normalize(bool is_full = true);

  // Replace uniformly-refined patch with single node
  void optimize(bool is_full = true);

  // Create an image of levels
  float * create_image (int n);

  // Return the number of levels
  int levels() { return levels_; }

 private:

  int levels_;
  Node4 * root_;

};
