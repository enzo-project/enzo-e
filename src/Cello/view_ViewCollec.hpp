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
#include <utility> // std::in_place_type_t
#include <variant>
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

private: // attributes

  std::variant<detail::SingleAddressViewCollec_<T>,
               detail::ArrOfPtrsViewCollec_<T>> collec_;

public: /// public interface
  /// default constructor
  ViewCollec() = default;

  /// construct a container of arrays without wrapping existing arrays
  ViewCollec(std::size_t n_arrays, const std::array<int,3>& shape) noexcept
    : collec_(std::in_place_type_t<detail::SingleAddressViewCollec_<T>>(),
              n_arrays, shape)
  { }

  /// construct from vector of arrays
  ViewCollec(const std::vector<CelloView<T, 3>>& v) noexcept
    : collec_(std::in_place_type_t<detail::ArrOfPtrsViewCollec_<T>>(), v)
  { }

  /// Destructor
  ~ViewCollec() = default;

  /// copy/move constructors & assignment operations
  ViewCollec(const ViewCollec& other)  = default;
  ViewCollec& operator=(const ViewCollec &other) = default;
  ViewCollec (ViewCollec&& other) = default;
  ViewCollec& operator=(ViewCollec&& other) = default;

  /// swaps the contents of `a` with `b`
  friend void swap(ViewCollec<T>& a, ViewCollec<T>& b) noexcept{
    a.collec_.swap(b.collec_);
  }

  /// Returns the number of arrays held in the collection
  std::size_t size() const noexcept {
    return std::visit([](auto&& arg) { return arg.size(); }, collec_);
  }

  /// Returns a shallow copy of the specified array
  const CelloView<T,3> operator[](std::size_t i) const noexcept {
    std::size_t n_elements = size();
    ASSERT2("ViewCollec::operator[]",
            "Can't retrieve value at %zu when size() is %zu",
            i, n_elements, i < n_elements);
    return std::visit([=](auto&& arg) { return arg[i]; }, collec_);
  }

  /// Indicates if the contained arrays are stored in a single 4D array
  /// (Alternatively they can be stored as an array of pointers)
  bool contiguous_items() const noexcept {
    return std::visit([](auto&& arg) { return arg.contiguous_items(); },
                      collec_);
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
    return std::visit([](auto&& arg) { return arg.get_backing_array(); },
                      collec_);
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
    return (*this)[0].shape(dim);
  }

  /// Return a new collection of subsections of each view in the collection
  ///
  /// @param slc_z, slc_y, slc_x Instance of CSlice that specify that are
  ///    passed to the `subarray` method of each contained array.
  ViewCollec<T> subarray_collec(const CSlice &slc_z,
                                const CSlice &slc_y,
                                const CSlice &slc_x) const noexcept
  {
    ViewCollec<T> out;
    auto fn = [&out, &slc_z, &slc_y, &slc_x](auto&& arg)
      { out.collec_ = arg.subarray_collec(slc_z,slc_y,slc_x); };
    std::visit(fn, collec_);
    return out;
  }

};

#endif /* VIEW_C_ARR_COLLEC_HPP */
