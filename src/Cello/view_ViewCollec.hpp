// See LICENSE_CELLO file for license and copyright information

/// @file     view_ViewCollec.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri Nov 26 2021
/// @brief    Declaration and implementation of the ViewCollec class template
//
// This class is primarily intended to be used to help implement the ArrayMap

#include <array>
#include <limits>
#include <memory> // std::shared_ptr, std::default_delete
#include <utility> // std::swap
#include <vector>

#ifndef VIEW_C_ARR_COLLEC_HPP
#define VIEW_C_ARR_COLLEC_HPP

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
                                    const std::vector<CelloView<T,3>> &arrays){
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
  struct SharedBuffer_{

    SharedBuffer_() = default;

    SharedBuffer_(const std::vector<T> &v) noexcept
      : arr_(std::shared_ptr<T>(new T[v.size()](), std::default_delete<T[]>())),
        length_(v.size())
    { for (std::size_t i = 0; i < length_; i++) { arr_.get()[i] = v[i]; } }

    const T& operator[](std::size_t i) const noexcept { return arr_.get()[i]; }
    T& operator[](std::size_t i) noexcept { return arr_.get()[i]; }
    std::size_t size() const noexcept {return length_;}

    void swap(SharedBuffer_<T>& other) noexcept
    {
      this->arr_.swap(other.arr_);
      std::swap(this->length_, other.length_);
    }

  private:
    std::shared_ptr<T> arr_;
    std::size_t length_;
  };

//----------------------------------------------------------------------

  template<typename T>
  class ArrOfPtrsViewCollec_{
    /// equivalent to an array of pointers (where each pointer is a CelloView)
    /// - we track the arrays lifetime with a std::shared_ptr to makes copies
    ///   cheaper. The disadvantage to this should be minimal since we don't
    ///   allow the "array" to be directly mutated after construction (note:
    ///   the "pointers" can't be mutated but the elements at the "pointer"
    ///   addresses can be mutated)

  public:

    ArrOfPtrsViewCollec_() = default;

    ArrOfPtrsViewCollec_(const std::vector<CelloView<T,3>> &arrays) noexcept
      : arrays_(arrays)
    { confirm_shared_shape_("ArrOfPtrsViewCollec_", arrays); }

    const CelloView<T,3> operator[](std::size_t index) const noexcept
    { return arrays_[index]; }

    std::size_t size() const noexcept { return arrays_.size(); }

    constexpr bool contiguous_items() const {return false;}

    const CelloView<T, 4> get_backing_array() const noexcept{
      ERROR("ArrOfPtrsViewCollec_::get_backing_array",
            "This is an invalid method call");
    }

    friend void swap(ArrOfPtrsViewCollec_<T>& a, ArrOfPtrsViewCollec_<T>& b)
      noexcept { a.arrays_.swap(b.arrays_); }

    // There might be some benefit to not directly constructing subarrays and
    // lazily evaluating them later
    ArrOfPtrsViewCollec_<T> subarray_collec(const CSlice &slc_z,
                                            const CSlice &slc_y,
                                            const CSlice &slc_x) const noexcept
    {
      std::vector<CelloView<T,3>> temp;
      for (std::size_t i = 0; i < arrays_.size(); i++){
        temp.push_back(arrays_[i].subarray(slc_z, slc_y, slc_x));
      }
      return ArrOfPtrsViewCollec_<T>(temp);
    }

  private:
    // ordered list of arrays
    SharedBuffer_<CelloView<T,3>> arrays_;
  };

//----------------------------------------------------------------------

  template<typename T>
  class SingleAddressViewCollec_{
    /// Implements an ordered collection of arrays with a single 4D CelloView.
    /// (the location of all of the contents are specified by a single address)

  public:
    SingleAddressViewCollec_() = default;

    SingleAddressViewCollec_(const std::size_t n_arrays,
                             const std::array<int,3>& shape) noexcept
      : backing_array_(int_coerce_(n_arrays), shape[0], shape[1], shape[2])
    { }
  
    const CelloView<T,3> operator[](std::size_t index) const noexcept
    { return backing_array_.subarray(int_coerce_(index)); }

    std::size_t size() const noexcept
    { return (backing_array_.is_null()) ? 0 : backing_array_.shape(0); }

    constexpr bool contiguous_items() const {return true;}

    const CelloView<T, 4> get_backing_array() const noexcept
    { return backing_array_; }

    friend void swap(SingleAddressViewCollec_<T>& a,
                     SingleAddressViewCollec_<T>& b) noexcept
    { swap(a.backing_array_, b.backing_array_); }

    SingleAddressViewCollec_<T> subarray_collec(const CSlice &slc_z,
                                                const CSlice &slc_y,
                                                const CSlice &slc_x) const
      noexcept
    {
      SingleAddressViewCollec_<T> out;
      out.backing_array_ = backing_array_.subarray(CSlice(0, nullptr),
                                                   slc_z, slc_y, slc_x);
      return out;
    }

  private:
    /// this holds the individual array elements
    CelloView<T, 4> backing_array_;
  };
}

template<typename T>
class ViewCollec{
  /// @class    ViewCollec
  /// @ingroup  View
  /// @brief    [\ref View] represents a collection of CelloViews of a
  ///           constant shape. While elements of these arrays can be mutated,
  ///           the arrays themselves can't be overwritten
  ///
  /// This primarily exists to help implement an ArrayMap
  ///
  /// This is a hybrid data-type that abstracts 2 separate implementations for
  /// collections of CelloViews. Under the hood, the arrays are either stored
  /// in a single 4D Contiguous Array (see SingleAddressViewCollec_) or an
  /// array of 3D CelloViews (see ArrOfPtrsViewCollec_)
  ///
  /// The implementation of this hybrid data-structure is fairly crude. If we
  /// decide that we want to support this hybrid data-type in the long-term,
  /// we probably want to consider using boost::variant/a backport of
  /// std::variant or take advantage of the fact that we can represent the
  /// contents of a SingleAddressViewCollec_ with a ArrOfPtrsViewCollec_

private:
  // tags how the arrays are stored (whether they're stored in a single
  // contiguous CelloView or stored like an array of pointers)
  enum class Tag {CONTIG, ARR_OF_PTR};
  Tag tag_;

  union{
    detail::SingleAddressViewCollec_<T> single_arr_;
    detail::ArrOfPtrsViewCollec_<T> arr_of_ptr_;
  };

private: // helper methods:

  /// default construct the active member
  void activate_member_(Tag tag, bool currently_initialized = true){
    if (currently_initialized) { dealloc_active_member_(); }
    tag_ = tag;

    switch (tag_){
      case Tag::CONTIG: {
        new(&single_arr_) detail::SingleAddressViewCollec_<T>;
        break;
      }
      case Tag::ARR_OF_PTR: {
        new(&arr_of_ptr_) detail::ArrOfPtrsViewCollec_<T>;
        break;
      }
    }
  }

  // used in destructor and to help switch the active member of the enum
  // need to explicitly call destructor because we use placement new
  void dealloc_active_member_() noexcept{
    // note, we need to use the full qualified names of destructors to avoid
    // errors with the nvidia compiler
    switch (tag_){
      case Tag::CONTIG: {
        single_arr_.detail::SingleAddressViewCollec_<T>::~SingleAddressViewCollec_();
        break;
      }
      case Tag::ARR_OF_PTR: {
        arr_of_ptr_.detail::ArrOfPtrsViewCollec_<T>::~ArrOfPtrsViewCollec_();
        break;
      }
    }
  }

public: /// public interface
  /// default constructor
  ViewCollec() noexcept { activate_member_(Tag::CONTIG, false); }

  /// construct a container of arrays without wrapping existing arrays
  ViewCollec(std::size_t n_arrays, const std::array<int,3>& shape) noexcept {
    activate_member_(Tag::CONTIG, false);
    single_arr_ = detail::SingleAddressViewCollec_<T>(n_arrays, shape);
  }

  /// construct from vector of arrays
  ViewCollec(const std::vector<CelloView<T, 3>>& v) noexcept{
    activate_member_(Tag::ARR_OF_PTR, false);
    arr_of_ptr_ = detail::ArrOfPtrsViewCollec_<T>(v);
  }

  /// copy constructor
  ViewCollec(const ViewCollec& other) noexcept{
    activate_member_(other.tag_, false);
    switch (other.tag_){
      case Tag::CONTIG:     { single_arr_ = other.single_arr_; break; }
      case Tag::ARR_OF_PTR: { arr_of_ptr_ = other.arr_of_ptr_; break; }
    }
  }

  /// copy assignment
  ViewCollec& operator=(const ViewCollec &other) noexcept {
    ViewCollec tmp(other);
    swap(*this,tmp);
    return *this;
  }

  /// move constructor
  ViewCollec (ViewCollec&& other) noexcept : ViewCollec(){ swap(*this,other); }

  /// move assignment
  ViewCollec& operator=(ViewCollec&& other) noexcept {
    swap(*this, other);
    return *this;
  }

  /// Destructor
  ~ViewCollec(){dealloc_active_member_();}

  /// swaps the contents of `a` with `b`
  friend void swap(ViewCollec<T>& a, ViewCollec<T>& b) noexcept{
    if (&a == &b) { return; }
    if (a.tag_ == b.tag_){
      switch (a.tag_){
        case Tag::CONTIG:     { swap(a.single_arr_, b.single_arr_); break; }
        case Tag::ARR_OF_PTR: { swap(a.arr_of_ptr_, b.arr_of_ptr_); break; }
      }
    } else if (a.tag_ == Tag::CONTIG){
      detail::SingleAddressViewCollec_<T> temp_single_arr;
      swap(temp_single_arr, a.single_arr_);
      a.activate_member_(Tag::ARR_OF_PTR, true);
      swap(a.arr_of_ptr_, b.arr_of_ptr_);
      b.activate_member_(Tag::CONTIG, true);
      swap(temp_single_arr, b.single_arr_);
    } else {
      swap(b,a);
    }
  }

  /// Returns the number of arrays held in the collection
  std::size_t size() const noexcept {
    switch (tag_){
      case Tag::CONTIG:     { return single_arr_.size(); }
      case Tag::ARR_OF_PTR: { return arr_of_ptr_.size(); }
    }
    ERROR("ViewCollec::size", "problem with exhaustive switch-statement");
  }

  /// Returns a shallow copy of the specified array
  const CelloView<T,3> operator[](std::size_t i) const noexcept {
    std::size_t n_elements = size();
    ASSERT2("ViewCollec::operator[]",
            "Can't retrieve value at %zu when size() is %zu",
            i, n_elements, i < n_elements);
    switch (tag_){
      case Tag::CONTIG:     { return single_arr_[i]; }
      case Tag::ARR_OF_PTR: { return arr_of_ptr_[i]; }
    }
    ERROR("ViewCollec::operator[]", "problem with exhaustive switch-statement");
  }

  /// Indicates if the contained arrays are stored in a single 4D array
  /// (Alternatively they can be stored as an array of pointers)
  bool contiguous_items() const noexcept {
    switch (tag_){
      case Tag::CONTIG:     { return single_arr_.contiguous_items(); }
      case Tag::ARR_OF_PTR: { return arr_of_ptr_.contiguous_items(); }
    }
    ERROR("ViewCollec::contiguous_items",
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
  const CelloView<T,4> get_backing_array() const noexcept {
    switch (tag_){
      case Tag::CONTIG:     { return single_arr_.get_backing_array(); }
      case Tag::ARR_OF_PTR: { return arr_of_ptr_.get_backing_array(); }
    }
    ERROR("ViewCollec::get_backing_array",
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
      ERROR("ViewCollec::array_shape", "Invalid when size() == 0");
    } else if (dim > 3) {
      ERROR("ViewCollec::array_shape", "invalid dimension");
    }
    return (*this)[0].shape(dim); // could be optimized based on tag
  }

  /// Return a new collection of subsections of each array in the collection
  ///
  /// @param slc_z, slc_y, slc_x Instance of CSlice that specify that are
  ///    passed to the `subarray` method of each contained array.
  ViewCollec<T> subarray_collec(const CSlice &slc_z,
                                const CSlice &slc_y,
                                const CSlice &slc_x) const noexcept
  {

    // initialize out and activate the union member matching this->tag_.
    ViewCollec<T> out;
    out.activate_member_(tag_, true);

    switch (tag_){
      case Tag::CONTIG: {
        out.single_arr_ = single_arr_.subarray_collec(slc_z,slc_y,slc_x);
        return out;
      }
      case Tag::ARR_OF_PTR: {
        out.arr_of_ptr_ = arr_of_ptr_.subarray_collec(slc_z,slc_y,slc_x);
        return out;
      }
    }
    ERROR("ViewCollec::subarray_collec_",
          "problem with exhaustive switch-statement");
  }

};

#endif /* VIEW_C_ARR_COLLEC_HPP */
