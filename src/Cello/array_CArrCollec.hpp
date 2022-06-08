// See LICENSE_CELLO file for license and copyright information

/// @file     array_CArrCollec.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri Nov 26 2021
/// @brief    Declaration and implementation of the CelloArray class template
//
// This class is primarily intended to be used to help implement the ArrayMap

#include <array>
#include <limits>
#include <vector>

#ifndef ARRAY_C_ARR_COLLEC_HPP
#define ARRAY_C_ARR_COLLEC_HPP

namespace detail {

  inline int int_coerce_(std::size_t arg) noexcept{
    if (arg > static_cast<std::size_t>(std::numeric_limits<int>::max())){
      ERROR1("coerce_to_int_", "Can't convert %zu to int", arg);
    }
    return static_cast<int>(arg);
  }

//----------------------------------------------------------------------

  template<typename T>
  inline void confirm_shared_shape_(const char* func_name,
                                    const std::vector<CelloArray<T,3>> &arrays){
    if (arrays.size() > 0){
      int ref_mz = arrays[0].shape(0);
      int ref_my = arrays[0].shape(1);
      int ref_mx = arrays[0].shape(2);
      for (std::size_t i = 1; i < arrays.size(); i++){
        int mz = arrays[i].shape(0);
        int my = arrays[i].shape(1);
        int mx = arrays[i].shape(2);

        ASSERT8(func_name,
                ("The shapes, (mz, my, mx), of arrays[%zu] and arrays[%zu] are "
                 "(%d, %d, %d) and (%d, %d, %d). The shapes must be the same"),
                0, i,
                ref_mz, ref_my, ref_mx, mz, my, mx,
                ((ref_mz == mz) && (ref_my == my) && (ref_mx == mx)));
      }
    }
  }

//----------------------------------------------------------------------

  template<typename T>
  class VecBackedCArrCollec_{
    /// equivalent to an array of pointers

  public:

    VecBackedCArrCollec_() = default;

    VecBackedCArrCollec_(const std::vector<CelloArray<T,3>> &arrays) noexcept
      : arrays_(arrays)
    { confirm_shared_shape_("VecBackedCArrCollec_",arrays); }

    const CelloArray<T,3> operator[](std::size_t index) const noexcept
    { return arrays_[index]; }

    std::size_t size() const noexcept { return arrays_.size(); }

    constexpr bool contiguous_items() const {return false;}

    const CelloArray<T, 4> get_backing_array() const noexcept{
      ERROR("VecBackedCArrCollec__::get_backing_array",
            "This is an invalid method call");
    }

    friend void swap(VecBackedCArrCollec_<T>& a, VecBackedCArrCollec_<T>& b)
      noexcept { a.arrays_.swap(b.arrays_); }

    // There might be some benefit to not directly constructing subarrays and
    // lazily evaluating them later
    VecBackedCArrCollec_<T> subarray_collec(const CSlice &slc_z,
                                            const CSlice &slc_y,
                                            const CSlice &slc_x) const noexcept
    {
      std::vector<CelloArray<T,3>> temp;
      for (std::size_t i = 0; i < arrays_.size(); i++){
        temp.push_back(arrays_[i].subarray(slc_z, slc_y, slc_x));
      }
      return VecBackedCArrCollec_<T>(temp);
    }

  private:
    // ordered list of arrays
    std::vector<CelloArray<T,3>> arrays_;
  };

//----------------------------------------------------------------------

  template<typename T>
  class SingleAddressCArrCollec_{
    /// Implements an ordered collection of arrays with a single 4D CelloArray.
    /// (the location of all of the contents are specified by a single address)

  public:
    SingleAddressCArrCollec_() = default;

    SingleAddressCArrCollec_(const std::size_t n_arrays,
                             const std::array<int,3>& shape) noexcept
      : backing_array_(int_coerce_(n_arrays), shape[0], shape[1], shape[2])
    { }
  
    const CelloArray<T,3> operator[](std::size_t index) const noexcept
    { return backing_array_.subarray(int_coerce_(index)); }

    std::size_t size() const noexcept
    { return (backing_array_.is_null()) ? 0 : backing_array_.shape(0); }

    constexpr bool contiguous_items() const {return true;}

    const CelloArray<T, 4> get_backing_array() const noexcept
    { return backing_array_; }

    friend void swap(SingleAddressCArrCollec_<T>& a,
                     SingleAddressCArrCollec_<T>& b) noexcept
    { swap(a.backing_array_, b.backing_array_); }

    SingleAddressCArrCollec_<T> subarray_collec(const CSlice &slc_z,
                                                const CSlice &slc_y,
                                                const CSlice &slc_x) const
      noexcept
    {
      SingleAddressCArrCollec_<T> out;
      out.backing_array_ = backing_array_.subarray(CSlice(0, nullptr),
                                                   slc_z, slc_y, slc_x);
      return out;
    }

  private:
    /// this holds the individual array elements
    CelloArray<T, 4> backing_array_;
  };
}

template<typename T>
class CArrCollec{
  /// @class    CArrCollec
  /// @ingroup  Array
  /// @brief    [\ref Array] represents a collection of CelloArrays of a
  ///           constant shape. While elements of these arrays can be mutated,
  ///           the arrays themselves can't be overwritten
  ///
  /// This primarily exists to help implement an ArrayMap
  ///
  /// This is a hybrid data-type that abstracts 2 separate implementations for
  /// collections of CelloArrays. Under the hood, the arrays are either stored
  /// in a single 4D Contiguous Array (see SingleAddressCArrCollec_) or a
  /// vector of 3D arrays (see VecBackedCArrCollec_)
  ///
  /// The implementation of this hybrid data-structure is fairly crude. If we
  /// decide that we want to support this hybrid data-type in the long-term,
  /// we probably want to consider using boost::variant/a backport of
  /// std::variant or take advantage of the fact that we can represent the
  /// contents of a SingleAddressCArrCollec_ with a VecBackedCArrCollec_

private:
  // tags how the arrays are stored (whether they're stored in a single
  // contiguous CelloArray or a vector of CelloArrays)
  enum class Tag {CONTIG, VEC};
  Tag tag_;

  union{
    detail::SingleAddressCArrCollec_<T> single_arr_;
    detail::VecBackedCArrCollec_<T> vec_of_arr_;
  };

private: // helper methods:

  /// default construct the active member
  void activate_member_(Tag tag, bool currently_initialized = true){
    if (currently_initialized) { dealloc_active_member_(); }
    tag_ = tag;

    using namespace detail;
    switch (tag_){
      case Tag::CONTIG: {new(&single_arr_) SingleAddressCArrCollec_<T>; break; }
      case Tag::VEC:    {new(&vec_of_arr_) VecBackedCArrCollec_<T>;     break; }
    }
  }

  // used in destructor and to help switch the active member of the enum
  // need to explicitly call destructor because we use placement new
  void dealloc_active_member_() noexcept{
    switch (tag_){
      case Tag::CONTIG: { single_arr_.~SingleAddressCArrCollec_(); break; }
      case Tag::VEC:    { vec_of_arr_.~VecBackedCArrCollec_(); break; }
    }
  }

public: /// public interface
  /// default constructor
  CArrCollec() noexcept { activate_member_(Tag::CONTIG, false); }

  /// construct a container of arrays without wrapping existing arrays
  CArrCollec(std::size_t n_arrays, const std::array<int,3>& shape) noexcept {
    activate_member_(Tag::CONTIG, false);
    single_arr_ = detail::SingleAddressCArrCollec_<T>(n_arrays, shape);
  }

  /// construct from vector of arrays
  CArrCollec(const std::vector<CelloArray<T, 3>>& v) noexcept{
    activate_member_(Tag::VEC, false);
    vec_of_arr_ = detail::VecBackedCArrCollec_<T>(v);
  }

  /// copy constructor
  CArrCollec(const CArrCollec& other) noexcept{
    activate_member_(other.tag_, false);
    switch (other.tag_){
      case Tag::CONTIG: { single_arr_ = other.single_arr_; break; }
      case Tag::VEC:    { vec_of_arr_ = other.vec_of_arr_; break; }
    }
  }

  /// copy assignment
  CArrCollec& operator=(const CArrCollec &other) noexcept {
    CArrCollec tmp(other);
    swap(*this,tmp);
    return *this;
  }

  /// move constructor
  CArrCollec (CArrCollec&& other) noexcept : CArrCollec(){ swap(*this,other); }

  /// move assignment
  CArrCollec& operator=(CArrCollec&& other) noexcept {
    swap(*this, other);
    return *this;
  }

  /// Destructor
  ~CArrCollec(){dealloc_active_member_();}

  /// swaps the contents of `a` with `b`
  friend void swap(CArrCollec<T>& a, CArrCollec<T>& b) noexcept{
    if (&a == &b) { return; }
    if (a.tag_ == b.tag_){
      switch (a.tag_){
        case Tag::CONTIG: { swap(a.single_arr_, b.single_arr_); break; }
        case Tag::VEC:    { swap(a.vec_of_arr_, b.vec_of_arr_); break; }
      }
    } else if (a.tag_ == Tag::CONTIG){
      detail::SingleAddressCArrCollec_<T> temp_single_arr;
      swap(temp_single_arr, a.single_arr_);
      a.activate_member_(Tag::VEC, true);
      swap(a.vec_of_arr_, b.vec_of_arr_);
      b.activate_member_(Tag::CONTIG, true);
      swap(temp_single_arr, b.single_arr_);
    } else {
      swap(b,a);
    }
  }

  /// Returns the number of arrays held in the collection
  std::size_t size() const noexcept {
    switch (tag_){
      case Tag::CONTIG: { return single_arr_.size(); }
      case Tag::VEC:    { return vec_of_arr_.size(); }
    }
    ERROR("CArrCollec::size", "problem with exhaustive switch-statement");
  }

  /// Returns a shallow copy of the specified array
  const CelloArray<T,3> operator[](std::size_t i) const noexcept {
    std::size_t n_elements = size();
    ASSERT2("CArrCollec::operator[]",
            "Can't retrieve value at %zu when size() is %zu",
            i, n_elements, i < n_elements);
    switch (tag_){
      case Tag::CONTIG: { return single_arr_[i]; }
      case Tag::VEC:    { return vec_of_arr_[i]; }
    }
    ERROR("CArrCollec::operator[]", "problem with exhaustive switch-statement");
  }

  /// Indicates if the contained arrays are stored in a single 4D array
  /// (Alternatively they can be stored as an array of pointers)
  bool contiguous_items() const noexcept {
    switch (tag_){
      case Tag::CONTIG: { return single_arr_.contiguous_items(); }
      case Tag::VEC:    { return vec_of_arr_.contiguous_items(); }
    }
    ERROR("CArrCollec::contiguous_items",
          "problem with exhaustive switch-statement");
  }

  /// Returns a shallow copy of the 4D array holding each contained array.
  ///
  /// The `n`th 3D subarray of `collec.get_backing_array()` is always a perfect
  /// alias of `collec[n]` (for any non-negative `n` less than `collec.size()`)
  ///
  /// @note
  /// The program will abort if this method is called on an object for which
  /// the `contiguous_items()` method returns `false`.
  const CelloArray<T,4> get_backing_array() const noexcept {
    switch (tag_){
      case Tag::CONTIG: { return single_arr_.get_backing_array(); }
      case Tag::VEC:    { return vec_of_arr_.get_backing_array(); }
    }
    ERROR("CArrCollec::get_backing_array",
          "problem with exhaustive switch-statement");
  }

  /// Returns the length along a given dimension of each contained array
  ///
  /// @param dim Indicates the dimension for which we want the shape. This
  ///    should be 0, 1, or 2
  ///
  /// @note
  /// The program will fail if this method is invoked and the collection holds
  /// zero elements.
  int array_shape(unsigned int dim) const noexcept
  {
    if (size() == 0) {
      ERROR("CArrCollec::array_shape", "Invalid when size() == 0");
    } else if (dim > 3) {
      ERROR("CArrCollec::array_shape", "invalid dimension");
    }
    return (*this)[0].shape(dim); // could be optimized based on tag
  }

  /// Return a new collection of subsections of each array in the collection
  ///
  /// @param slc_z, slc_y, slc_x Instance of CSlice that specify that are
  ///    passed to the `subarray` method of each contained array.
  CArrCollec<T> subarray_collec(const CSlice &slc_z,
                                const CSlice &slc_y,
                                const CSlice &slc_x) const noexcept
  {

    // initialize out and activate the union member matching this->tag_.
    CArrCollec<T> out;
    out.activate_member_(tag_, true);

    switch (tag_){
      case Tag::CONTIG: {
        out.single_arr_ = single_arr_.subarray_collec(slc_z,slc_y,slc_x);
        return out;
      }
      case Tag::VEC:    {
        out.vec_of_arr_ = vec_of_arr_.subarray_collec(slc_z,slc_y,slc_x);
        return out;
      }
    }
    ERROR("CArrCollec::subarray_collec_",
          "problem with exhaustive switch-statement");
  }

};

#endif /* ARRAY_C_ARR_COLLEC_HPP */
