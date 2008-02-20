#ifndef CLASS_ARRAY
#define CLASS_ARRAY 
// $Id$
/**
 * @file
 * @brief Header file for the Array class
 * @author James Bordner 
 * @version 1.0
 *
 * Attributes
 *
 * (*) n_  int[4]
 * (*) N_  int
 * (*) a_  Scalar * 
 *
 * Operations
 *
 *     public:
 *
 * (*) Array()
 * (*) Array(int  n0, int  n1=1, int  n2=1, int n3=1)
 * (*) ~Array()
 * (*) void copy (Array *)
 * (*) void resize (int   n0, int   n1=1, int   n2=1, int   n3=1)
 * (*) void size (int * n0, int * n1=0, int * n2=0, int * n3=0) const
 * (*) int  length() const
 * (*) Scalar * values () const
 * ( ) Scalar & operator() (int  i0, int  i1=1, int  i2=1, int i3=1)
 *
 *     private:
 *
 * (*) void allocate_ (int n0, int   n1=1, int   n2=1, int   n3=1))
 * (*) void deallocate_()
 * (*) void copy_ (Scalar * a)
 *
 */
// $Log$
 
//----------------------------------------------------------------------

/// Create a new Array object
 
/**
 */
 
class Array {

  // PRIVATE ATTRIBUTES

private:

  /// Length of array: N_ == n_[0]*n_[1]*n_[2]*n_[3]

  int     N_;

  /// Shape of array, right-padded with 1's

  int     n_[4];

  /// Array values stored in column-major ordering

  Scalar *a_;

  // PUBLIC OPERATIONS

public:

  Array();
  Array(int  n0, int  n1=1, int  n2=1, int n3=1);
  ~Array();
  void copy (const Array &);
  void resize (int   n0, int   n1=1, int   n2=1, int   n3=1);
  void size (int * n0, int * n1=0, int * n2=0, int * n3=0) const;
  int  length() const;
  Scalar * values () const;
  Scalar & operator() (int  i0, int  i1=0, int  i2=0, int i3=0);

  // PRIVATE OPERATIONS

private:

  /// Allocate array values
  void allocate_ (int n0, int n1=1, int n2=1, int n3=1);

  /// Deallocate array values
  void deallocate_();

  /// Copy array values a[0:N-1] to Array
  void copy_ (Scalar * a);

};

#endif
