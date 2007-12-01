//345678901234567890123456789012345678901234567890123456789012345678901234567890

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

class Grid;
typedef Grid * pGrid;
//typedef int[12] pEntries;


// Maximum number of nonstencil nonzero elements in any row
// limiting case: coarse zone in corner bounded 
// by fine grid patches on three faces

const int MAX_NONZERO = 12; 

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
    _adjacent_covered_,
    _last_ =  _adjacent_covered_
  };

  static const int _num_types_;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Label[] ENTRIES SHOULD MATCH LabelName VALUES IN faces.cpp
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  static const char * LabelName[];

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

 private:

  // data defined at grid creation

  /// Arrays for face-zone types

  Label * label_[3][2];

  /// Arrays for zone neighbors

  pGrid * neighbor_[3][2];

  /// Counters for nonstencil element indices required by
  /// HYPRE_SStructMatrixAddToValues() in "entries".  Initialized
  /// to be the number of stencil entries.

  int * entry_coarse_[3][2];    
  int * entry_fine_[3][2];    

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

  /// Set or get the face-zone type of the i,jth zone on the given
  /// axis (0 to 2) and face (0 to 1).  No error checking is performed
  /// on axis, face, i or j.

  Label &label (int axis, int face, int i, int j) throw ()
  { return label_[axis][face][i+n1_[axis]*j]; };


  /// Set the face-zone type of all zones on the given axis (0
  /// to 2) and face (0 to 1).  No error checking is performed on
  /// axis or face.

  void label (int axis, int face, Label label) throw ()
  { 
    for (int i=0; i<n1_[axis]; i++) {
      for (int j=0; j<n2_[axis]; j++) {
	label_[axis][face][i+n1_[axis]*j] = label; 
      }
    }
  };

  /// Set or get the facing grid neighbor of the i,jth zone on the given
  /// axis (0 to 2) and face (0 to 1).  No error checking is performed
  /// on axis, face, i or j.

  pGrid & neighbor (int axis, int face, int i, int j) throw ()
  { return neighbor_[axis][face][i+n1_[axis]*j]; };


  /// Set the facing grid neighbor of all zones on the given axis (0
  /// to 2) and face (0 to 1).  No error checking is performed on
  /// axis or face.

  void neighbor (int axis, int face, pGrid neighbor) throw ()
  { 
    for (int i=0; i<n1_[axis]; i++) {
      for (int j=0; j<n2_[axis]; j++) {
	neighbor_[axis][face][i+n1_[axis]*j] = neighbor; 
      }
    }
  };

  /// Set or get the current entry index of the facing grid neighbor
  /// of the i,jth zone on the given axis (0 to 2) and face (0 to 1).
  /// No error checking is performed on axis, face, i or j.

  int & entry_fine (int axis, int face, int i, int j) throw ();
  int & entry_coarse (int axis, int face, int i, int j) throw ();

  /// Set the facing grid neighbor of all zones on the given axis (0
  /// to 2) and face (0 to 1).  No error checking is performed on
  /// axis or face.

  void entry_coarse (int axis, int face, int value) throw ()
  { 
    for (int i=0; i<n1_[axis]; i++) {
      for (int j=0; j<n2_[axis]; j++) {
	entry_coarse_[axis][face][i+n1_[axis]*j] = value; 
      }
    }
  };
  void entry_fine (int axis, int face, int value) throw ()
  { 
    for (int i=0; i<n1_[axis]; i++) {
      for (int j=0; j<n2_[axis]; j++) {
	entry_fine_[axis][face][i+n1_[axis]*j] = value; 
      }
    }
  };

  int n1 (int axis) { return n1_[axis]; };
  int n2 (int axis) { return n2_[axis]; };

  void print() throw();

private:

  void alloc_ (int *n) throw ();
  void dealloc_ () throw ();
  int & entry_ (int axis, int face, int i, int j, int * entry [3][2]) throw ();

  //--------------------------------------------------------------------

};
