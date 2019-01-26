
#ifndef ENZO_ENZO_ARRAY_HPP
#define ENZO_ENZO_ARRAY_HPP
#include <stdio.h>

// Planned changes:
//   - Behavior of operator() if fewer arguments are specified than the number
//     of dimensions. The current behavior is problematic if EnzoArray is a
//     view of a subarray of a different array. This will ultimately be
//     addressed by declaring the number of dimensions as a template argument
//     and then comparing the number of indices to the number of dimensions
//     (this should be done at compile time)
//   - Memory can leak for a view of an EnzoArray if the original owns the
//     underlying pointer and is deallocated. To address this we start tracking
//     views of other EnzoArrays that own pointers, and reorganizing ownership
//     of the pointer upon deallocation. This does not address lack of memory
//     safety for EnzoArrays wrapping arbitrary pointers (the only possible way
//     to fix that would be to restrict EnzoArrays to wrap pointers from Field
//     or a subclass and then have field notify the EnzoArray if a pointer is
//     deallocated - Going to avoid doing this).
//   - Probably going to make initialization and Construction occur in one
//     step, rather than 2 steps. To help facillitate this, move and copy
//     constructors will be implemented

// Implementation Notes:
//   - Main motivation: use overloaded operator() for multi-dimensional access
//     of C style arrays
//   - If this object is adopted for use on a larger scale, then it may be
//     better to define this object in the Cello Directory
//   - This should not be confused with the EnzoArray in the original enzo code.
//     It may be better to rename this CelloArray
//   - Supports up to 4D arrays (everything is implicitly a 4D array)
//   - Like Athena++ arrays (and unlike numpy style arrays) EnzoArrays will
//     always return a single value. Imagine a 3D array with shape {mz,my,mx}
//       - For both types of arrays: array(k,j,i) is equivalent to accessing 
//         index ((k*my + j)*mx + i) of a contiguous array
//       - For EnzoArrays: array(j,i) is equivalent to index (j*mx + i)
//                         array(i) is equivalent to index i
//         Both of these cases would return subarrays for numpy style arrays
//       - Users should avoid this behavior. If EnzoArray is a view of a
//         subarray, unexpected behavior will occur.
//   - dimensions are ordered by increasing indexing speed of each axis. A 4d
//     array has dimensions: (dim3, dim2, dim1, dim0)
//   - To allow for writing functions performing directional grid operations
//     that can be generalized to arbitrary grid direction, EnzoArray is
//     allowed to represent a subsection of another array. This is equivalent
//     to saying a[1:3,0:4,0:4] for a numpy array.
//       - For example, separate EnzoArrays may be used to represent the left
//         and right interface states that are centered on cell faces along an
//         arbitrary axis. In this case, the for-loop body could remain the
//         same and different subarrays could be operated on. 
//       - This means that the EnzoArray will have a slightly larger memory
//         footprint
//       - The alternative would be to permute the axes order of the array.
//         This is undesirable for 2 reasons:
//           1. Generalizing functions to arbitrary dimensions would be more
//              difficult. The simplest way to do so would involve iterating
//              over indices in different orders
//           2. A fast implementation would be more complicated and a secondary
//              class definition may be necessary for speed.


// Defining the macro called CHECK_BOUNDS macro means that the bounds of an
// EnzoArray, are checked everytime operator() is called or a subarray is
// initialized.

#ifdef CHECK_BOUNDS
#  define CHECK_BOUND1D(i,shape_)                                             \
  ASSERT("EnzoArray","Invalid index", i<shape_[0]);
#  define CHECK_BOUND2D(j,i,shape_)					      \
  ASSERT("EnzoArray", "Invalid index", i<shape_[0] && j<shape_[1]);
#  define CHECK_BOUND3D(k,j,i,shape_)                                         \
  ASSERT("EnzoArray", "Invalid index",                                        \
	 i<shape_[0] && j<shape_[1] && k<shape_[2]);
#  define CHECK_BOUND4D(l,k,j,i,shape_)					\
  ASSERT("EnzoArray", "Invalid index",                                        \
	 i<shape_[0] && j<shape_[1] && k<shape_[2] && l<shape_[3]);
#else
#  define CHECK_BOUND1D(i,shape_)       /* ... */
#  define CHECK_BOUND2D(j,i,shape_)     /* ... */
#  define CHECK_BOUND3D(k,j,i,shape_)   /* ... */
#  define CHECK_BOUND4D(l,k,j,i,shape_) /* ... */
#endif

template<typename T>
class EnzoArray
{

public:
  // Constructors of an EnzoArray by wrapping an existing c style array
  EnzoArray()
    : data_(NULL),
      offset_(0),
      shape_(),
      stride_(),
      owns_ptr_(false)
  { };

  // Destructor
  ~EnzoArray();

  // Initialize an EnzoArray and allocate memory
  void initialize(int dim3, int dim2, int dim1, int dim0);

  // Overload Initialize for fewer dimensions
  void initialize(int dim2, int dim1, int dim0){
    initialize(1,dim2,dim1,dim0);}
  void initialize(int dim1, int dim0){
    initialize(1,1,dim1,dim0);}
  void initialize(int dim0){
    initialize(1,1,1,dim0);}

  // Initialize EnzoArray by wrapping an existing c style array
  void initialize_wrapper(T* array, int dim3, int dim2, int dim1, int dim0);

  // Overload InitializeWrapper for fewer dimensions
  void initialize_wrapper(T* array, int dim2, int dim1, int dim0){
    initialize_wrapper(array, 1, dim2, dim1, dim0); }
  void initialize_wrapper(T* array, int dim1, int dim0){
    initialize_wrapper(array, 1, 1, dim1, dim0); }
  void initialize_wrapper(T* array, int dim0){
    initialize_wrapper(array, 1, 1, 1, dim0); }

  // Initialize EnzoArray that represents a subsection of an existing EnzoArray
  void initialize_subarray(EnzoArray<T> &array, int dim3_start, int dim3_stop,
			   int dim2_start, int dim2_stop, int dim1_start,
			   int dim1_stop, int dim0_start, int dim0_stop);

  void initialize_subarray(EnzoArray<T> &array, int dim2_start, int dim2_stop,
			   int dim1_start, int dim1_stop, int dim0_start,
			   int dim0_stop){
    initialize_subarray(array, 0, array.shape_[3], dim2_start, dim2_stop,
			dim1_start, dim1_stop, dim0_start, dim0_stop);}
  void initialize_subarray(EnzoArray<T> &array, int dim1_start, int dim1_stop,
			  int dim0_start, int dim0_stop){
    initialize_subarray(array, 0, array.shape_[3], 0, array.shape_[2],
			dim1_start, dim1_stop, dim0_start, dim0_stop);}
  void initialize_subarray(EnzoArray<T> &array, int dim0_start, int dim0_stop){
    initialize_subarray(array, 0, array.shape_[3], 0, array.shape_[2],
			0, array.shape_[1], dim0_start, dim0_stop);}

  // Access array Elements
  T &operator() (const int i){
    CHECK_BOUND1D(i,shape_);
    return data_[offset_+i]; }
  T operator() (const int i) const{
    CHECK_BOUND1D(i,shape_);
    return data_[offset_+i]; }

  T &operator() (const int j, const int i){
    CHECK_BOUND2D(j,i,shape_);
    return data_[offset_ + i + j*stride_[0]]; }
  T operator() (const int j, const int i) const{
    CHECK_BOUND2D(j,i,shape_);
    return data_[offset_ + i + j*stride_[0]];}

  T &operator() (const int k, const int j, const int i){
    CHECK_BOUND3D(k,j,i,shape_);
    return data_[offset_ + i + j*stride_[0] + k*stride_[1]]; }
  T operator() (const int k, const int j, const int i) const{
    CHECK_BOUND3D(k,j,i,shape_);
    return data_[offset_ + i + j*stride_[0] + k*stride_[1]]; }

  T &operator() (const int l, const int k, const int j, const int i){
    CHECK_BOUND4D(l,k,j,i,shape_);
    return data_[offset_ + i + j*stride_[0] + k*stride_[1] + l*stride_[2]]; }
  T operator() (const int l, const int k, const int j, const int i) const{
    CHECK_BOUND4D(l,k,j,i,shape_);
    return data_[offset_ + i + j*stride_[0] + k*stride_[1] + l*stride_[2]]; }

  int length_dim0() {return shape_[0];}
  int length_dim1() {return shape_[1];}
  int length_dim2() {return shape_[2];}
  int length_dim3() {return shape_[3];}

  int size(){
    int out = 1;
    for (int i=0; i<4; i++){
      out*=shape_[i];
    }
    return out;
  }

  // This function exists for debugging.
  const T* get_ptr(){return data_;}

protected: // methods
  void init_helper_(T* array, int dim3, int dim2, int dim1, int dim0,
		    int offset, bool owns_ptr);

  // Deallocates memory if necessary
  void cleanup_helper_();

private: // attributes

  // array data
  T *data_;

  // offset of the first element from the start of the pointer
  int offset_;

  // shape_ lists values as dim0, dim1,...
  // dim0 is the length of the axis of fastest indexing.
  // Indices are listed in reverse Order: (...,dim2,dim1,dim0)
  int shape_[4];

  // stride_ lists values as dim1, dim2,...
  // length of the stride for dim0 is always 1.
  int stride_[3];

  // Whether or not the object is repsonsible for freeing the pointer
  bool owns_ptr_;
};

template<typename T>
void EnzoArray<T>::cleanup_helper_()
{
  if (owns_ptr_){
    delete []data_;
  }
  owns_ptr_ = false;
}

template<typename T>
EnzoArray<T>::~EnzoArray()
{
  cleanup_helper_();
}

template<typename T>
void EnzoArray<T>::init_helper_(T* array, int dim3, int dim2, int dim1,
				int dim0, int offset, bool owns_ptr)
{
  // to prevent memory leaks, need to be sure memory is deallocated if an array
  // is allocated multiple times
  if (owns_ptr_){
    if (array != data_){
      cleanup_helper_();
    } else {
      // If this instance of EnzoArray already owned the pointer, it needs to
      // maintain ownership
      owns_ptr = true;
    }
  }

  data_ = array;
  owns_ptr_ = owns_ptr;
  offset_ = offset;

  shape_[0] = dim0; shape_[1] = dim1; shape_[2] = dim2; shape_[3] = dim3;
  stride_[0] = dim0;
  stride_[1] = dim1*stride_[0];
  stride_[2] = dim2*stride_[1];
}

// initialize by allocating data
template<typename T>
void EnzoArray<T>::initialize(int dim3, int dim2, int dim1, int dim0)
{
  T* array = new T[dim3*dim2*dim1*dim0](); // allocate and set entries to 0
  init_helper_(array, dim3, dim2, dim1, dim0, 0, true);
}

// initialize by wrapping existing c-style array
template<typename T>
void EnzoArray<T>::initialize_wrapper(T* array, int dim3, int dim2, int dim1,
				      int dim0)
{
  init_helper_(array, dim3, dim2, dim1, dim0, 0, false);
}

// initialize EnzoArray that represents a subsection of an existing EnzoArray
template<typename T>
void EnzoArray<T>::initialize_subarray(EnzoArray<T> &array,
				       int dim3_start, int dim3_stop,
				       int dim2_start, int dim2_stop,
				       int dim1_start, int dim1_stop,
				       int dim0_start, int dim0_stop)
{
  CHECK_BOUND4D(dim3_start,dim2_start,dim1_start,dim0_start,array.shape_);
  CHECK_BOUND4D(dim3_stop-1,dim2_stop-1,dim1_stop-1,dim0_stop-1,array.shape_);
  
  int offset = (array.offset_ + dim0_start +
		dim1_start*array.stride_[0] +
		dim2_start*array.stride_[1] +
		dim3_start*array.stride_[2]);
  
  int dim3 = dim3_stop - dim3_start;
  int dim2 = dim2_stop - dim2_start;
  int dim1 = dim1_stop - dim1_start;
  int dim0 = dim0_stop - dim0_start;

  init_helper_(array.data_, dim3, dim2, dim1, dim0, offset, false);
  stride_[0] = array.stride_[0];
  stride_[1] = array.stride_[1];
  stride_[2] = array.stride_[2];
}

#endif /* ENZO_ENZO_ARRAY_HPP */
