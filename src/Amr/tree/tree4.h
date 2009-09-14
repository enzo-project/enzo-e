#define MAX_LEVELS 80

class Node4;

class Tree4 {

 public:

  // create a tree refined to the array non-zeros
  // assume array width and height are powers of 2
  Tree4();

  // Refine down to array
  void refine(const bool * array, int nd0, int nd1);

  // print levels
  void print_levels();

  // Refine down to array
  void normalize();

  // Create an image of levels
  float * create_image (int * nd0, int * nd1, int n0, int n1, int cell_width);

  // Return the number of levels
  int levels() { return levels_; }
 private:

  int levels_;
  Node4 * root_;

};
