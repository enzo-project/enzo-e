

class Node4;

class Tree4 {

 public:

  // create a tree refined to the array non-zeros
  // assume array width and height are powers of 2
  Tree4(const bool *array, int width, int height);

  /// Print the tree
  void print(); 


 private:

  Node4 * root_;

  const bool * array_;
  int width_;
  int height_;
};
