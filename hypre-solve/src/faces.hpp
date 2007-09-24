
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

public:

  /// Types of face-zones

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Label[] ENTRIES SHOULD MATCH LabelName VALUES IN faces.cpp
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  enum Label {
    _first_,
    _unknown_ = _first_,
    _boundary_,
    _coarse_,
    _fine_,
    _neighbor_,
    _covered_,
    _last_ =  _covered_
  };

  static const char * LabelName[];

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

 private:

  // data defined at grid creation

  /// Boolean arrays for faces; values are relative level offsets

  Label * label_[3][2];

  /// Leading dimension of arrays
  int  n1_[3]; 
  int  n2_[3]; 
  int  n_[3]; 

  //--------------------------------------------------------------------

 public:

  /// Categorization of face zones

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

  /// Relative offset of neighboring cell outside the grid on the
  /// given axis and face.  No error checking on axis, face, i or j.

  Label &label (int axis, int face, int i, int j) throw ()
  { return label_[axis][face][i+n1_[axis]*j]; };

  void label (int axis, int face, Label label) throw ()
  { 
    for (int i1=0; i1<n1_[axis]; i1++) {
      for (int i2=0; i2<n2_[axis]; i2++) {
	label_[axis][face][i1+n1_[axis]*i2] = label; 
      }
    }
  };

  int n1 (int axis) { return n1_[axis]; };
  int n2 (int axis) { return n2_[axis]; };

  void print() throw();

private:

  void alloc_ (int *n) throw ();
  void dealloc_ () throw ();

  //--------------------------------------------------------------------

};


