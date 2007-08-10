
/// Discret class include file

/**
 * 
 * The Discret class stores and operates on discretization aspects of a grid.
 * 
 *   Each Grid class has a corresponding Discret class.  The Discret
 *   class is used to set up the matrix elements for the unified 
 *   AMR grid.
 *  
 *   Setting up the matrix elements is done using the following
 *   steps:
 *
 *   1. Define the stencil for all interior grid elements
 *
 *   2. Zero-out matrix elements covered by a refined grid
 *
 *   3. Handle matrix elements connecting neighboring grids
 *   - Handle matrix elements defining the boundary conditions
 *
 *   - Handle coarse unknowns adjacent to fine unknowns in child
 *   - Handle fine unknowns adjacent to coarse unknowns in parent
 *
 *   - Handle coarse unknowns adjacent to fine unknowns in neighbor's child
 *   - Handle fine unknowns adjacent to coarse unknowns in parent's neighbor
 *     
 *   - Handle any remaining connections
 *
 *   Note that the matrix generated is not symmetric.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

class Grid;

class Discret
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

public:
  enum enum_neighbor_cell {_unknown_,_coarse_, _fine_, _same_, _boundary_} ;

 private:

  // data defined at grid creation

  /// Boolean arrays for faces
  enum_neighbor_cell * neighbor_cell_[3][2];
  /// Leading dimension of arrays
  int  n1_[3]; 
  int  n2_[3]; 
  int  n_[3]; 

  //--------------------------------------------------------------------

 public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  /// Create a discretization for a grid
  Discret (int *n) throw ();

  ~Discret () throw ();

  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  // Data access

  /// Type of cell (_unknown_,_coarse_,_fine_,or _same_) away from grid on given axis/face
  enum_neighbor_cell &neighbor_cell (int axis, int face, int i, int j) throw ()
  { return neighbor_cell_[axis][face][i+n1_[axis]*j]; };

  void print() throw();

private:

  void alloc_ (int *n) throw ();
  void dealloc_ () throw ();

  //--------------------------------------------------------------------

};


