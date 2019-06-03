// See LICENSE_CELLO file for license and copyright information

/// @file     array_CelloArray.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 30 2019
/// @brief    Declaration and implementation of the CelloArray class template

#include <stdio.h>
#include <cstddef>
#include <type_traits>
#include <limits>
#include <memory>

#ifndef ARRAY_CELLO_ARRAY_HPP
#define ARRAY_CELLO_ARRAY_HPP

//----------------------------------------------------------------------

// CelloArray
// Like arrays in Athena++, the indices are listed in order of increasing
//   access speed. Imagine a 3D array with shape {mz,my,mx} array(k,j,i) is
//   equivalent to accessing index ((k*my + j)*mx + i) of the pointer
// Dimensions are numbered with increasing indexing speed (dim0, dim1, ...)

//----------------------------------------------------------------------
/// @typedef intp
/// @brief  Array indexing type
///
/// This convention is borrowed from numpy. We set it the larger precision 
/// integer type: int OR ptrdiff_t to guarantee int promotion
/// (Alternatively, we could define it as int)

typedef std::conditional<sizeof(int) >= sizeof(std::ptrdiff_t),
			   int, std::ptrdiff_t>::type intp;

//----------------------------------------------------------------------
/// @def    ARRAY_SIZE_MAX
/// @brief  Maximum number of array elements
///
/// If we forced intp to always be int, then we'd use min(INT_MAX,PTRDIFF_MAX)

#define ARRAY_SIZE_MAX PTRDIFF_MAX

//----------------------------------------------------------------------


class CSlice
{
  /// @class    CSlice
  /// @ingroup  Array
  /// @brief    [\ref Array] represents a slice along a single axis of a
  ///           CelloArray
  ///
  /// To specify a specify a slice including indices from start through stop-1
  /// call CSlice(start, stop). Negative indexing is also supported.
  ///
  /// To indicate that the slice should extend to the end of the axis, use:
  /// - CSlice(start, NULL) OR CSlice(start, nullptr)
  /// To indicate that the slice should extend from the start of the axis to
  /// some known stop index, use:
  /// - CSlice(0, stop) OR CSlice(NULL, stop) OR CSlice(nullptr, stop)
  /// To indicate that the slice should extend over the full axis, use:
  /// - CSlice(0, NULL), CSlice(NULL, NULL), CSlice(0, nullptr) OR
  ///   CSlice(0, nullptr)
  ///
  /// 2 important usage notes:
  /// 1. Due to the implementation defined nature of NULL, 0 can be used in
  ///    place of NULL or nullptr. However the convention is to avoid this
  /// 2. The class is equipped with a default constructor to allow for
  ///    construction of arrays of CSlice. However, all default constructed
  ///    instances MUST be assigned the value of non-default constructed
  ///    instances before passing them to CelloArray.subarray (This is done to
  ///    help catch mistakes)

public:
  
  // ints should always be implicitly promoted
  CSlice(const intp start, const intp stop)
    : start_(start), stop_(stop), initialized_(true)
  {
    if ((start >= 0 && stop > 0) || (start < 0 && stop < 0)){
      ASSERT("CSlice", "start must be less than stop.", stop>start);
    }
  }

  // Exists only to allow construction of arrays of slices. The constructed
  // instance must be assigned a non-default constructed value before use
  CSlice() : start_(0), stop_(0), initialized_(false) {}

  /*
  CSlice(const intp start, std::nullptr_t stop)
    : start_(start), stop_(0), initialized_(true)
  { }

  CSlice(std::nullptr_t start, const intp stop)
    : start_(0), stop_(stop), initialized_(true)
  { }

  CSlice(std::nullptr_t start, std::nullptr_t stop)
    : start_(0), stop_(0), initialized_(true)
  { }
  */

  // Gets the start index of the slice
  intp get_start() const
  {
    ASSERT("CSlice", ("Default constructed CSlices cannot be used without "
		      "explicit assignment of values."), initialized_);
    return start_;
  }

  /// Returns the stop index of the slice. If the stop index should be the
  /// length of sliced axis, 0 is returned.
  intp get_stop() const
  {
    ASSERT("CSlice", ("Default constructed CSlices cannot be used without "
		      "explicit assignment of values."), initialized_);
    return stop_;
  }
    
private:
  /// start index
  intp start_;
  /// stop index (Should be 0 if the full length if the dimension will be used)
  intp stop_;

  /// true if constructed with non-default constructor OR if assigned values
  /// from a non-default constructed instance.
  bool initialized_;
};


// Defining the macro called CHECK_BOUNDS macro means that the bounds of an
// EnzoArray, are checked everytime operator() is called
template<typename T>
bool check_bounds_(std::size_t *shape, T first) {return *shape > first;}
template<typename T, typename... Rest>
bool check_bounds_(std::size_t *shape, T first, Rest... rest){
  //(get's unrolled at compile time)
  return (*shape > first) && check_bounds_(++shape, rest...);
}

#ifdef CHECK_BOUNDS
#  define CHECK_BOUND3D(shape, k, j, i)                                       \
  ASSERT("FixedDimArray_","Invalid index", check_bounds_(shape,k,j,i));
#  define CHECK_BOUNDND(shape, ARGS)                                          \
  ASSERT("FixedDimArray_","Invalid index", check_bounds_(shape, ARGS...));
#else
#  define CHECK_BOUND3D(shape, k, j, i)   /* ... */
#  define CHECK_BOUNDND(shape, ARGS)      /* ... */
#endif


// Defining the macro called CHECK_FINITE_ELEMENTS returns an error everytime a
// NaN or inf gets retrieved from the array or an array is initialized by
// wrapping an existing pointer holding a NaN or inf
template<typename T>
bool check_if_finite_(const T elem){
  return (std::is_floating_point<T>::value) ? std::isfinite(elem) : true;
}

#ifdef CHECK_FINITE_ELEMENTS
#  define CHECK_IF_FINITE(value)                                              \
  ASSERT("FixedDimArray_","Non-Finite Value", check_if_finite_(value));
#  define CHECK_IF_ARRAY_FINITE(value)    this->assert_all_entries_finite_();
#else
#  define CHECK_IF_FINITE(value)          /* ... */
#  define CHECK_IF_ARRAY_FINITE(value)    /* ... */
#endif



// To define CelloArray to have arbitrary dimensions it needs to accept a
// variable number of arguments to indicate shape during construction, to
// produce a subarray, and to access elements. The number of dimensions of the
// array as a template argument and accept values with variadic template
// arguments to check that the appropriate the number of values are specified
// at compile time. In each case, we need to guaruntee that all arguments
// are a given type. The solution is based on:
//   - https://stackoverflow.com/a/28253503
//   - https://stackoverflow.com/a/31767710
template <bool...> struct bool_pack;
template <bool... vals> using all_true = std::is_same<bool_pack<true, vals...>,
						      bool_pack<vals..., true>>;
// May want to add more possible int types (e.g. long,short,etc.)
#define REQUIRE_TYPE(T,type1)                                                \
  class = std::enable_if<all_true<(std::is_same<T,type1>::value)...>::value>
#define REQUIRE_TYPE2(T,type1,type2)					     \
  class = std::enable_if<all_true<(std::is_same<T,type1>::value ||           \
				   std::is_same<T,type2>::value)...>::value>
#define REQUIRE_INT(T) REQUIRE_TYPE2(T,intp,int)


// Helper function to compute pointer index
// T is expected to be intp or int
template<typename T>
intp calc_index_(intp* stride, T first){return first;}

template<typename T>
intp calc_index_(intp* stride, T first, T last){
  return (*stride)*first + last;
}

template<typename T, typename... Rest>
intp calc_index_(intp* stride, T first, Rest... rest){
  // gets unrolled at compile time
  return (*stride)*first + calc_index_(++stride, rest...);
}

// The following 2 functions are helpful while dynamically iterating
template<typename T>
intp calc_index_(const std::size_t D, const intp offset,
		 const intp* stride, const T* indices){
  std::size_t out = offset + indices[D-1];
  for (std::size_t i = 0; i+1<D;i++){
    out += stride[i]*indices[i];
  }
  return out;
}

// outer means indices for "outer" loops
inline void increment_outer_indices_(std::size_t D, intp *indices, intp *shape,
				     bool &continue_outer_iter){
  std::size_t i = D-1;
  while (0 != (i--)){
    indices[i]+=1;
    if (indices[i] != shape[i]){
      return;
    } else if (i > 0){
      indices[i] = 0;
    }
  }
  continue_outer_iter = false;
}


template<typename T>
class dataWrapper
{
  // tracks the underlying CelloArray pointer and ownership of it
public:

  dataWrapper(T* data, bool owns_ptr) : data_(data), owns_ptr_(owns_ptr) { };
  ~dataWrapper() {if (owns_ptr_) { delete[] data_; }}
  T* get() const noexcept { return data_; }

private:
  T *data_;
  bool owns_ptr_;
};


template<typename T, std::size_t D>
class CelloArray;

template<typename T, std::size_t D>
class TempArray_;

template<typename T, std::size_t D>
class FixedDimArray_
{
public:

  // Destructor
  ~FixedDimArray_() { cleanup_helper_();}

  // operator() - used to access array Elements
  template<typename... Args, REQUIRE_INT(Args)>
  T &operator() (Args... args) {
    static_assert(D==sizeof...(args),
		  "Number of indices don't match number of dimensions");
    CHECK_BOUNDND(shape, args)
    CHECK_IF_FINITE(data_[offset_ + calc_index_(stride_,args...)]);
    return data_[offset_ + calc_index_(stride_,args...)];
  }
  template<typename... Args, REQUIRE_INT(Args)>
  T operator() (Args... args) const {
    static_assert(D==sizeof...(args),
		  "Number of indices don't match number of dimensions");
    CHECK_BOUNDND(shape, args)
    CHECK_IF_FINITE(data_[offset_ + calc_index_(stride_,args...)]);
    return data_[offset_ + calc_index_(stride_,args...)];
  }

  // Specialized implementation for 3D arrays (reduces compile time)
  T &operator() (const int k, const int j, const int i){
    static_assert(D==3, "3 indices should only be specified for 3D arrays");
    CHECK_BOUND3D(shape, k, j, i)
    CHECK_IF_FINITE(data_[offset_ + k*stride_[0] + j*stride_[1] + i])
    return data_[offset_ + k*stride_[0] + j*stride_[1] + i];
  }
  T operator() (const int k, const int j, const int i) const{
    static_assert(D==3, "3 indices should only be specified for 3D arrays");
    CHECK_BOUND3D(shape, k, j, i)
    CHECK_IF_FINITE(data_[offset_ + k*stride_[0] + j*stride_[1] + i])
    return data_[offset_ + k*stride_[0] + j*stride_[1] + i];
  }


  // Produce a copy of the array.
  TempArray_<T,D> deepcopy();
  
  // Returns a subarray with same D. Expects an instance of CSlice for each
  // dimension. For a 3D CelloArray, the function declaration might look like:
  // CelloArray<T,3> subarray(CSlice k_slice, CSlice j_slice, CSlice i_slice);
  //
  // The Special Case of no arguments returns the full array
  template<typename... Args, REQUIRE_TYPE(Args,CSlice)>
  TempArray_<T,D> subarray(Args... args);

  int shape(unsigned int dim){
    ASSERT1("FixedDimArray_", "%ui is greater than the number of dimensions",
	    dim, dim<D);
    return (int)shape_[dim];
  }

  intp size(){
    intp out = 1;
    for (std::size_t i=0; i<D; i++){
      out*=shape_[i];
    }
    return out;
  }

  // Only arrays with the same numbers of dimensions can be swapped
  friend void swap(FixedDimArray_<T,D> &first, FixedDimArray_<T,D> &second){
    std::swap(first.dataMgr_, second.dataMgr_);
    std::swap(first.data_, second.data_);
    std::swap(first.offset_, second.offset_);
    std::swap(first.shape_, second.shape_);
    std::swap(first.stride_, second.stride_);
  }

protected: // methods to be reused by subclasses

  FixedDimArray_()
    : dataMgr_(),
      data_(NULL),
      offset_(0),
      shape_(),
      stride_()
  { };

  // Construct a numeric array that allocates its own data
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  FixedDimArray_(Args... args);

  // Construct a numeric array that wraps an existing pointer
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  FixedDimArray_(T* array, Args... args);
  
  void init_helper_(std::shared_ptr<dataWrapper<T>> &dataMgr, intp shape[D],
		    intp offset){
    data_ = dataMgr->get();
    dataMgr_ = dataMgr;
    offset_ = offset;

    std::size_t i = D;
    while (i>0){
      --i;
      shape_[i] = shape[i];
      if (i + 1 == D){
	stride_[i] = 1;
      } else {
	stride_[i] = shape_[i+1] * stride_[i+1];
      }
    }
  }

  void cleanup_helper_(){ data_ = NULL; }

  // this method is for debugging only, it is not made available to users
  void assert_all_entries_finite_();

public: // attributes
  // (these are only public so that various subclasses can access these)
  std::shared_ptr<dataWrapper<T>> dataMgr_; // manages ownership of data_
  // pointer to data (copied from dataMgr to provide faster access to elements)
  T* data_;
  intp offset_; // offset of the first element from the start of the pointer
  intp shape_[D]; // lists dimensions with increasing indexing speed
  intp stride_[D]; // stride_[D-1] is always 1
};


// check that the array shape is allowed (all positive and not too big)
inline void check_array_shape_(intp shape[], std::size_t D)
{
  ASSERT("FixedDimArray_", "Positive dimensions are required.", shape[0]>0);
  ASSERT1("FixedDimArray_", "The array cannot exceed %ld elements.",
	  (long)ARRAY_SIZE_MAX, (ARRAY_SIZE_MAX / shape[0]) < 1);
  intp cur_size = 1;

  for (std::size_t i = 0; i+1 < D; i++){
    cur_size *= shape[i];
    ASSERT("FixedDimArray_", "Positive dimensions are required.",shape[i+1]>0);
    ASSERT1("FixedDimArray_", "The array cannot exceed %ld elements.",
	    (long)ARRAY_SIZE_MAX, (ARRAY_SIZE_MAX / shape[i+1]) < cur_size);
  }
}
  


// Constructor of CelloArray by allocating new data
template<typename T, std::size_t D>
template<typename... Args, class>
FixedDimArray_<T,D>::FixedDimArray_(Args... args)
{
  static_assert(D==sizeof...(args), "Incorrect number of dimensions");
  intp shape[D] = {((intp)args)...};
  check_array_shape_(shape, D);

  intp size = 1;
  for (std::size_t i=0; i < D; i++){
    size *= shape[i];
  }
  T* data = new T[size](); // allocate and set entries to 0

  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(data,true);
  init_helper_(dataMgr, shape, 0);
}


// Constructor of array that wraps an existing c-style array
template<typename T, std::size_t D>
template<typename... Args, class>
FixedDimArray_<T,D>::FixedDimArray_(T* array, Args... args)
{
  static_assert(D==sizeof...(args), "Incorrect number of dimensions");
  intp shape[D] = {((intp)args)...};
  check_array_shape_(shape, D);

  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(array,false);
  init_helper_(dataMgr, shape, 0);
  CHECK_IF_ARRAY_FINITE(value)
}


// Prepares an array of slices that refer to absolute start and stop values
// along each dimension. Also checks that the slices are valid
inline void prep_slices_(const CSlice* slices, const intp shape[],
			 const std::size_t D, CSlice* out_slices)
{
   for (std::size_t i=0; i<D; i++){
     intp start = slices[i].get_start();
     if (start < 0){
       start += shape[i]; // start was negative
     }

     intp stop = slices[i].get_stop();
     if (stop <= 0){
       // includes negative values of stop and case when it's equal to zero
       // (which means that the slice should stop at shape[i])
       stop += shape[i];
     }

     ASSERT3("FixedDimArray_",
	     "slice start of %ld doesn't lie in bound of dim %ld of size %ld.",
	     (long)slices[i].get_start(), (long)i, (long)shape[i],
	     start < shape[i]);
     ASSERT3("FixedDimArray_",
	     "slice stop of %d doesn't lie in bound of dim %ld of size %ld.",
	     (long)slices[i].get_stop(), (long)i, (long)shape[i],
	     stop <= shape[i]);
     ASSERT4("FixedDimArray_", ("slice stop (%ld) doesn't exceed slice start "
				"(%ld) must for dim %ld of size %ld."),
	     (long)slices[i].get_start(), (long)slices[i].get_stop(), (long)i,
	     (long)shape[i], stop>start);
     out_slices[i] = CSlice(start,stop);
  }
}

// Returnd TempArray_ representing a view of a subarray of the current instance
template<typename T, std::size_t D>
template<typename... Args, class>
TempArray_<T,D> FixedDimArray_<T,D>::subarray(Args... args){
  static_assert(D == sizeof...(args) || 0 == sizeof...(args),
		"Number of slices don't match number of dimensions");
  TempArray_<T,D> subarray;
  if (sizeof...(args) == 0) {
    subarray.init_helper_(dataMgr_,shape_,offset_);
  } else {
    CSlice in_slices[] = {args...};
    CSlice slices[D];
    prep_slices_(in_slices, shape_, D, slices);

    intp new_shape[D];
    intp new_offset = offset_;
    for (std::size_t dim=0; dim<D; dim++){
      new_shape[dim] = slices[dim].get_stop() - slices[dim].get_start();
      new_offset += slices[dim].get_start() * stride_[dim];
    }

    subarray.init_helper_(dataMgr_,new_shape,new_offset);
    for (std::size_t dim=0; dim<D; dim++){
      subarray.stride_[dim] = stride_[dim];
    }
  }
  return subarray;
}



template<typename T, std::size_t D>
TempArray_<T,D> FixedDimArray_<T,D>::deepcopy()
{
  T* data = new T[size()]; // allocate but don't initialize values
  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(data,true);
  TempArray_<T,D> out;
  out.init_helper_(dataMgr, shape_, 0);
  out = *this;
  return out;
}



template<typename T, std::size_t D>
void FixedDimArray_<T,D>::assert_all_entries_finite_()
{
  bool continue_outer_iter = true;
  intp indices[D] = {}; // all elements to 0
  while (continue_outer_iter){
    intp index = calc_index_(D, this->offset_, this->stride_, indices);
    for (intp i = 0; i<this->shape_[D-1]; i++){

      if (!check_if_finite_(this->data_[index])){
	//  let's get the indices as a string
	std::string str_indices = "";
	for (std::size_t j=0; j+1<D; j++){
	  str_indices+= std::to_string(indices[j]) + ", ";
	}
	str_indices += std::to_string(i);

	ASSERT1("FixedDimArray_", "The element at (%s) has a non_NaN value.",
		str_indices.c_str(),false);

      }
      index++;
    }
    increment_outer_indices_(D, indices, this->shape_,continue_outer_iter);
  }
}


template<typename T, std::size_t D>
class TempArray_ : public FixedDimArray_<T,D>
{
  friend class FixedDimArray_<T,D>;
  friend class CelloArray<T,D>;

public:
  // Assigns to each element of *this the value of val
  TempArray_<T,D>& operator=(const T& val);
  
  // Assigns to each element of *this the value of the corresponding element in
  // other. Sizes must match
  TempArray_<T,D>& operator=(const TempArray_<T,D>& other){
    assign_helper_(other.offset_,other.stride_,other.shape_, other.data_);
    return *this;
  }

  TempArray_<T,D>& operator=(const CelloArray<T,D>& other){
    assign_helper_(other.offset_,other.stride_, other.shape_, other.data_);
    return *this;
  }

private:
  TempArray_() : FixedDimArray_<T,D>() { }

  void assign_helper_(const intp o_offset, const intp *o_stride,
		      const intp *o_shape, const T* o_data){
    for (std::size_t i = 0; i<D; i++){
      ASSERT("TempArray_","shapes aren't the same.",this->shape_[i]==o_shape[i]);
    }
    bool continue_outer_iter = true;
    intp indices[D] = {}; // all elements to 0
    while (continue_outer_iter){
      intp index = calc_index_(D, this->offset_, this->stride_, indices);
      intp o_index = calc_index_(D, o_offset, o_stride, indices);
      for (intp i = 0; i<this->shape_[D-1]; i++){
	this->data_[index] = o_data[o_index];
	index++; o_index++;
      }
      increment_outer_indices_(D, indices, this->shape_, continue_outer_iter);
    }
  }

};


template<typename T, std::size_t D>
TempArray_<T,D>& TempArray_<T,D>::operator=(const T& val)
{
  bool continue_outer_iter = true;
  intp indices[D] = {}; // all elements to 0
  while (continue_outer_iter){
    intp index = calc_index_(D, this->offset_, this->stride_, indices);
    for (intp i = 0; i<this->shape_[D-1]; i++){
      this->data_[index] = val;
      index++;
    }
    increment_outer_indices_(D, indices, this->shape_,continue_outer_iter);
  }
  return *this;
}


template<typename T, std::size_t D>
class CelloArray : public FixedDimArray_<T,D>
{
public:
  // Default constructor. Constructs an unallocated CelloArray.
  CelloArray() : FixedDimArray_<T,D>() { }

  // Construct a numeric array that allocates its own data
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  CelloArray(Args... args) : FixedDimArray_<T,D>(args...) { }

  // Construct a numeric array that wraps an existing pointer
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  CelloArray(T* array, Args... args) : FixedDimArray_<T,D>(array, args...) { }

  // Copy constructor. Constructs a shallow copy of other.
  CelloArray(const CelloArray<T,D>& other){
    this->init_helper_(other.dataMgr_, other.shape_, other.offset_);
  }

  // Move constructor. Constructs the array with the contents of other
  CelloArray(CelloArray<T,D>&& other) : CelloArray() {swap(*this,other);}
  CelloArray(TempArray_<T,D>&& other) : CelloArray() {swap(*this,other);}

  // Copy assignment operator. Makes *this a shallow copy of other. (Contents
  // of shallow copies and subarrays of *this are unaffected)
  CelloArray<T,D>& operator=(const CelloArray<T,D>& other){
    this->cleanup_helper_();
    init_helper_(other.dataMgr_, other.shape_, other.offset_);
    return *this;
  }

  // Move assignment operator. Replaces the contents of *this with those of
  // other. (Contents of shallow copies and subarrays of *this are unaffected)
  CelloArray<T,D>& operator=(CelloArray<T,D>&& other) {
    swap(*this,other);
    return *this;
  }
  CelloArray<T,D>& operator=(TempArray_<T,D>&& other) {
    swap(*this,other);
    return *this;
  }
};

#endif /* ARRAY_CELLO_ARRAY_HPP */
