
/// Faces class include file

/**
 * 
 * The Faces class stores and operates on discretization aspects of a grid.
 * 
 * Each Grid class has a corresponding Faces class.  The Faces class
 * is used to help set up the matrix elements for the unified AMR
 * grid.
 *  
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

class Faces
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

 private:

  // data defined at grid creation

  /// Boolean arrays for faces; values are relative level offsets

  const static int _unknown_;
  const static int _boundary_;

  int * neighbor_cell_[3][2];

  /// Leading dimension of arrays
  int  n1_[3]; 
  int  n2_[3]; 
  int  n_[3]; 

  //--------------------------------------------------------------------

 public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  /// Create a Faces clas for a grid
  Faces (int *n) throw ();

  ~Faces () throw ();

  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  // Data access

  /// Relative offset of neighboring cell outside the grid on the given axis and face

  int &neighbor_cell (int axis, int face, int i, int j) throw ()
  { return neighbor_cell_[axis][face][i+n1_[axis]*j]; };

  void print() throw();

private:

  void alloc_ (int *n) throw ();
  void dealloc_ () throw ();

  //--------------------------------------------------------------------

};


