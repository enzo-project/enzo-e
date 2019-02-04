
#ifndef ENZO_ENZO_ARRAY_HPP
#define ENZO_ENZO_ARRAY_HPP

#include <stdio.h>
#include <cstddef>
#include <list>
#include <type_traits>
#include <memory>

// Planned changes:
//   - adding EnzoArray<T,D>::copy() to produce deep copies
// Possible changes:
//   - making a special TempArray class that represents the view that is
//     returned each time. These objects are meant to be deleted immediately
//     after use (and should not be constructable through normal means). The
//     idea would be to overload the operator= operation to allow for setting
//     many elements equal to one value AND for copying values values between
//     arrays.

// Implementation Notes:
//   - Main motivation: use overloaded operator() for multi-dimensional access
//     of C style arrays
//   - If this object is adopted for use on a larger scale, then it may be
//     better to define this object in the Cello Directory
//   - This should not be confused with the EnzoArray in the original enzo code.
//     It may be better to rename this CelloArray.
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
inline void check_bound_ND_(std::size_t N, int* indices, int* shape){
  // indices lists indices in order of decreasing dimension:
  //   (idim2, idim1, idim0)
  // shape_ list values in order of increasing dimensions:
  //   (dim0, dim1, dim2)
  for (std::size_t i=0; i<N; i++){
    ASSERT("EnzoArray", "Invalid index",  indices[i]<shape[N-i-1]);
  }
}

#ifdef CHECK_BOUNDS
#  define CHECK_BOUND1D(i,shape_)                                             \
  ASSERT("EnzoArray","Invalid index", i<shape_[0]);
#  define CHECK_BOUND2D(j,i,shape_)					      \
  ASSERT("EnzoArray", "Invalid index", i<shape_[0] && j<shape_[1]);
#  define CHECK_BOUND3D(k,j,i,shape_)                                         \
  ASSERT("EnzoArray", "Invalid index",                                        \
	 i<shape_[0] && j<shape_[1] && k<shape_[2]);
#  define CHECK_BOUNDND(N,ind,shape_) check_bound_ND_(N, ind, shape_);
  
#else
#  define CHECK_BOUND1D(i,shape_)       /* ... */
#  define CHECK_BOUND2D(j,i,shape_)     /* ... */
#  define CHECK_BOUND3D(k,j,i,shape_)   /* ... */
#  define CHECK_BOUNDND(N,ind,shape_)   /* ... */
#endif

// To define EnzoArray to have arbitrary dimension we need to accept a variable
// number of arguments to indicate shape during construction, to produce a
// subarray, and to access elements. In each case, it would be optimal to
// check that the appropriate the number of variables are specified at compile
// time. Consequently, we specify the number of dimensions of the array as a
// template argument, and in each case will use variadic template arguments.
// 
// We need to guaruntee that in each of the above functions that the arguments
// are all integers (or other appropriate types for indexing). The solution is
// based on:
//   - https://stackoverflow.com/a/28253503
//   - https://stackoverflow.com/a/31767710

// Use REQUIRE_INT_DECL(T) when declaring the class method.
template <bool...> struct bool_pack;
template <bool... vals> using all_true = std::is_same<bool_pack<true, vals...>,
						      bool_pack<vals..., true>>;
// May want to add more possible int types (e.g. long,short,etc.)
#define REQUIRE_INT(T)                                                     \
  class = std::enable_if<all_true<(std::is_same<T,std::size_t>::value ||   \
				   std::is_same<T,int>::value)...>::value>





// dataWrapper tracks the underlying EnzoArray pointer and tracks ownership
// Since this is written with C++11 in mind - can't make a shared_ptr directly
// represent an array pointer

template<typename T>
class dataWrapper
{
public:

  dataWrapper(T* data, bool owns_ptr) : data_(data), owns_ptr_(owns_ptr) { };

  ~dataWrapper(){
    if (owns_ptr_){
      delete[] data_;
    }
  }

  T* get() const noexcept { return data_; }

private:
  T *data_;
  bool owns_ptr_;
};


// Here we use dataWrapper to wrap the pointer
// D is the number of dimensions of EnzoArray

template<typename T, std::size_t D>
class EnzoArray
{
public:
  // Default constructor. Constructs an unallocated EnzoArray.
  EnzoArray()
    : dataMgr_(),
      data_(NULL),
      offset_(0),
      shape_(),
      stride_()
  { };

  // Construct a numeric array that allocates its own data
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  EnzoArray(Args... args);

  // Construct a numeric array that wraps an existing pointer
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  EnzoArray(T* array, Args... args);

  // Copy constructor. Constructs a shallow copy of other.
  EnzoArray(const EnzoArray<T,D>& other){
    init_helper_(other.dataMgr_, other.shape_, other.offset_);
  }

  // Move constructor. Constructs the array with the contents of other using
  // move semantics
  EnzoArray(EnzoArray<T,D>&& other) : EnzoArray() {swap(*this,other);}

  // Destructor
  ~EnzoArray() { cleanup_helper_();}

  // Copy assignment operator. Makes *this a shallow copy of other. Data
  // previously owned by *this is forfeited before making the copy (the
  // contents of shallow copies and subarrays of *this are unaffected)
  EnzoArray<T,D>& operator=(const EnzoArray<T,D>& other){
    cleanup_helper_();
    init_helper_(other.dataMgr_, other.shape_, other.offset_);
    return *this;
  }

  // Move assignment operator. Replaces the contents of *this with those of
  // other. Data previously owned by *this is forfeited before making the copy
  // the contents of shallow copies and subarrays of *this are unaffected)
  EnzoArray<T,D>& operator=(EnzoArray<T,D>&& other) {
    swap(*this,other);
    return *this;
  }

  // operator() - used to access array Elements
  // General implementation

  template<typename... Args, REQUIRE_INT(Args)>
  T &operator() (Args... args) {
    return data_[calc_index_(args...)];}
  template<typename... Args, REQUIRE_INT(Args)>
  T operator() (Args... args) const {
    return data_[calc_index_(args...)];}


  // Specialized implementation for 1D, 2D, and 3D arrays
  T &operator() (const int i){
    static_assert(D==1, "1 index should only be specified for 1D arrays");
    CHECK_BOUND1D(i,shape_);
    return data_[offset_+i]; }
  T operator() (const int i) const{
    static_assert(D==1, "1 index should only be specified for 1D arrays");
    CHECK_BOUND1D(i,shape_);
    return data_[offset_+i]; }

  T &operator() (const int j, const int i){
    static_assert(D==2, "2 indices should only be specified for 2D arrays");
    CHECK_BOUND2D(j,i,shape_);
    return data_[offset_ + i + j*stride_[0]]; }
  T operator() (const int j, const int i) const{
    static_assert(D==2, "2 indices should only be specified for 2D arrays");
    CHECK_BOUND2D(j,i,shape_);
    return data_[offset_ + i + j*stride_[0]];}

  T &operator() (const int k, const int j, const int i){
    static_assert(D==3, "3 indices should only be specified for 3D arrays");
    CHECK_BOUND3D(k,j,i,shape_);
    return data_[offset_ + i + j*stride_[0] + k*stride_[1]]; }
  T operator() (const int k, const int j, const int i) const{
    static_assert(D==3, "3 indices should only be specified for 3D arrays");
    CHECK_BOUND3D(k,j,i,shape_);
    return data_[offset_ + i + j*stride_[0] + k*stride_[1]]; }


  // Return a subarray of the same dimension. Expects twice as many arguments
  // as there are dimensions (specifying the start and stop values along each
  // dimension)
  // 
  // For a 3D EnzoArray, the function declaration would look like:
  // EnzoArray<T,3> subarray(int dim2_start, int dim2_stop,
  //                         int dim1_start, int dim1_stop,
  //                         int dim0_start, int dim0_stop);
  template<typename... Args, REQUIRE_INT(Args)>
  EnzoArray<T, D> subarray(Args... args);

  int dim_size(unsigned int num){
    ASSERT("EnzoArray",
	   "get_dim should for dimension numbers less than dim.",
	   num<D);
    return shape_[num];
  }

  // Depreciated - the following needs to be removed
  int length_dim0() {return shape_[0];}
  int length_dim1() {return shape_[1];}
  int length_dim2() {return shape_[2];}
  int length_dim3() {return shape_[3];}

  int size(){
    int out = 1;
    for (int i=0; i<D; i++){
      out*=shape_[i];
    }
    return out;
  }

  // I believe this means that only arrays with the same numbers of dimensions
  // can be swapped
  friend void swap(EnzoArray<T,D> &first, EnzoArray<T,D> &second){
    std::swap(first.dataMgr_, second.dataMgr_);
    std::swap(first.data_, second.data_);
    std::swap(first.offset_, second.offset_);
    std::swap(first.shape_, second.shape_);
    std::swap(first.stride_, second.stride_);
  }

private: // methods
  // Helps initialize function
  void init_helper_(std::shared_ptr<dataWrapper<T>> &dataMgr, int shape[D],
		    int offset);

  // Deallocates memory if necessary
  void cleanup_helper_();

  // helps compute the index of underlying pointer
  template<typename... Args>
  int calc_index_(Args... args);

protected: // attributes

  // manages ownership of data_ pointer
  std::shared_ptr<dataWrapper<T>> dataMgr_;

  // pointer to array data
  T* data_;

  // offset of the first element from the start of the pointer
  int offset_;

  // shape_ lists values as dim0, dim1,...
  // dim0 is the length of the axis of fastest indexing.
  // Indices are listed in reverse Order: (...,dim2,dim1,dim0)
  int shape_[D];

  // stride_ lists values as dim1, dim2,...
  // length of the stride for dim0 is always 1.
  int stride_[D-1];
};


template<typename T, std::size_t D>
void EnzoArray<T,D>::cleanup_helper_()
{
  data_ = NULL;
}

template<typename T, std::size_t D>
void EnzoArray<T,D>::init_helper_(std::shared_ptr<dataWrapper<T>> &dataMgr,
				  int shape[D], int offset)
{
  data_ = dataMgr->get();
  dataMgr_ = dataMgr;
  offset_ = offset;

  for (std::size_t i=0; i<D;i++){
    shape_[i] = shape[i];
    if (i == 1){
      stride_[0] = shape_[0];
    } else if (i > 1){
      stride_[i-1] = shape_[i-1]*stride_[i-2];
    }
  }
}


// Constructor of EnzoArray by allocating new data
template<typename T, std::size_t D>
template<typename... Args, class>
EnzoArray<T,D>::EnzoArray(Args... args)
{
  static_assert(D==sizeof...(args),
		"The specified shape has a different number of dimensions "
		"than the template value.");
  int temp_shape[D] = {((int)args)...};
  int shape[D];
  for (std::size_t i=0; i < D; i++){
    ASSERT("EnzoArray", "The array shape must consis of positive elements.",
	   temp_shape[D-i-1]>0);
    shape[i] = temp_shape[D-i-1];
  }

  int size = 1;
  for (std::size_t i=0; i < D; i++){
    size *= shape[i];
  }
  
  T* data = new T[size](); // allocate and set entries to 0

  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(data,true);
  init_helper_(dataMgr, shape, 0);
}


// Constructor of EnzoArray that wraps an existing c-style array
template<typename T, std::size_t D>
template<typename... Args, class>
EnzoArray<T,D>::EnzoArray(T* array, Args... args){

  static_assert(D==sizeof...(args),
		"The specified shape has a different number of dimensions "
		"than the template value.");
  int temp_shape[D] = {((int)args)...};
  int shape[D];
  for (std::size_t i=0; i < D; i++){
    ASSERT("EnzoArray", "The array shape must consis of positive elements.",
	   temp_shape[D-i-1]>0);
    shape[i] = temp_shape[D-i-1];
  }

  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(array,false);
  init_helper_(dataMgr, shape, 0);
}

// Returnd an EnzoArray that represents a view of a subarray of the
// current instance
template<typename T, std::size_t D>
template<typename... Args, class>
EnzoArray<T,D> EnzoArray<T,D>::subarray(Args... args){
  static_assert(D+D == sizeof...(args),
		"subarray expects the number of arguments to equal twice the "
		"array's dimension.");
  int vals[2*D] = {((int)args)...};
  // vals lists values in order of decreasing dimension:
  //   (dim2_start,dim2_stop, dim1_start, dim1_stop, dim0_start, dim0_stop)
  // shape_ list values in order of increasing dimensions:
  //   (dim0, dim1, dim2)
  for (std::size_t i=0; i<D; i++){
    ASSERT2("EnzoArray<T,D>::subarray",
	    "The start value, %d, along dimension, %d, can't be negative.",
	    vals[2*i], (int)(D-1-i),vals[2*i]>=0);
    ASSERT3("EnzoArray<T,D>::subarray",
	    ("The start value, %d, must be less than the stop value, %d, along "
	     "dimension %d."), vals[2*i], vals[2*i+1], (int)(D-1-i),
	    vals[2*i]<vals[2*i+1]);
    ASSERT3("EnzoArray<T,D>::subarray",
	    "The stop value, %d, along dimension, %d, must be <= %d.",
	    vals[2*i+1], (int)(D-1-i),shape_[D-i-1],vals[2*i+1]<=shape_[D-i-1]);
  }


  int new_shape[D];
  int new_offset = offset_;
  for (std::size_t dim=0; dim<D; dim++){
    int shape_ind = (int)dim;
    int slice_start_ind = (int)(D-dim-1)*2;
    int slice_stop_ind = (int)(D-dim-1)*2+1;
    if (dim == 0){
      new_offset += vals[slice_start_ind];
    } else {
      new_offset += vals[slice_start_ind] * stride_[shape_ind-1];
    }
    new_shape[shape_ind] = vals[slice_stop_ind] - vals[slice_start_ind];
  }

  EnzoArray<T,D> subarray;

  subarray.init_helper_(dataMgr_,new_shape,new_offset);
  for (std::size_t dim=0; dim<D-1; dim++){
    subarray.stride_[dim] = stride_[dim];
  }
  return subarray;
}

// Helper function that computes the index of the underlying pointer
template<typename T, std::size_t D>
template<typename... Args>
int EnzoArray<T,D>::calc_index_(Args ...args)
{
  static_assert(D==sizeof...(args),
		"D indices are expected for a D-dimensional array.");
  int indices[D] = {((int)args)...};
  CHECK_BOUNDND(D,indices,shape_)
  int index = offset_ + indices[D-1]; 
  for (std::size_t dim = 1; dim<D; dim++){
    index += indices[D-1-dim] * stride_[dim-1];
  }
  return index;
}

#endif /* ENZO_ENZO_ARRAY_HPP */
