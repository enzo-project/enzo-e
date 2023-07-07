// See LICENSE_CELLO file for license and copyright information

/// @file     view_ViewMap.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon July 3 2023
/// @brief    Declaration and implementation of the ViewMap class template
///
/// If we ever decide that template instantiation is taking too much time, (at
/// the time of writing, this seems very unlikely after profiling compile-times
/// with clang's -ftime-trace), we can look into using extern template
///  - using extern template, we can force the template to be instantiated once
///    in a single translation for each floating-point and integer types
///  - while using extern template would be BAD for CelloView (since it
///    prevents inlining), there probably wouldn't be significant runtime cost
///    here.
///
/// All insertions occur in the constructor, and the properties of the
/// contained arrays can't be mutated. This choice:
///   - facillitates enforcement that all contained arrays have a fixed shape
///   - makes it easier to order the entries in an arbitrary order. This
///     facilitates optimizations in the Riemann Solver when the values are
///     initialized in the order expected by the Riemann Solver
///   - facillitates some optimizations that make instances relatively cheap to
///     copy.
///
/// If necessary, a number of optimizations could be made to the implementation
/// that might make key lookups faster. These optimizations could take
/// advantage of the following factors:
///    - Assumptions about the max key size and the max capacity of the map.
///      For example, if the max key size never exceeds ~22 characters and
///      there are never more than ~128 entries, it would probably be optimal
///      to store the strings in-place (improving cache locallity).
///    - It would also be worth considering whether linear search is faster
///      (since the arrays are small)

#ifndef VIEW_VIEW_MAP_HPP
#define VIEW_VIEW_MAP_HPP

#include <array>
#include <string>
#include <vector>

template<typename T>
class ViewMap {
  /// @class ViewMap
  /// @ingroup View
  /// @brief [\ref View] Stores instances of CelloView
  ///
  /// This was originally a concrete type called EnzoEFltArrayMap

public: // interface

  typedef T value_type;
  typedef typename std::add_const<T>::type const_value_type;
  typedef typename std::remove_const<T>::type nonconst_value_type;

  friend class ViewMap<const_value_type>;

  ViewMap(std::string name)
    : name_(name)
  { }

  ViewMap()
    : ViewMap("")
  { }

  /// Constructs a map that owns all of its own data.
  ViewMap(std::string name, const std::vector<std::string> &keys,
          const std::array<int,3>& shape);

  ViewMap(const std::vector<std::string> &keys,
          const std::array<int,3>& shape)
    : ViewMap("", keys, shape)
  { }

  /// Constructs a map that wraps existing data.
  ///
  /// Each array must have the same shape.
  ViewMap(std::string name, const std::vector<std::string> &keys,
          const std::vector<CelloView<T,3>> &arrays);

  ViewMap(const std::vector<std::string> &keys,
          const std::vector<CelloView<T,3>> &arrays)
    : ViewMap("", keys, arrays)
  { }

  /// conversion constructor that facilitates implicit casts from
  /// ViewMap<nonconst_value_type> to ViewMap<const_value_type>
  ///
  /// @note
  /// This is only defined for instances of ViewMap for which T is const-
  /// qualified. If it were defined in cases where T is not const-qualified,
  /// then it would duplicate the copy-constructor.
  template<class = std::enable_if<std::is_same<T, const_value_type>::value>>
  ViewMap(const ViewMap<nonconst_value_type> &other)
    : name_(other.name_),
      str_index_map_(other.str_index_map_),
      arrays_(other.arrays_)
  { }

  /// Returns a reference to the mapped array associated with the specified
  /// index/key
  ///
  /// Unlike the analogous ``std::map::operator[]`` method, this cannot be used
  /// to insert a new value. This behaves similarly to ``std::map::at``
  ///
  /// @note
  /// We explicitly return by value to prevent users from overwriting the
  /// CelloView instance stored internally (in other words the pointer to
  /// data). Such an operation would not work for one of the backends.
  ///
  /// @note
  /// Alternatively, we could return a constant reference, but that's like
  /// returning a constant reference to a pointer. In reality, that might be
  /// more advantageous here since a CelloView is larger than a generic ptr.
  CelloView<T, 3> operator[] (const std::string& key) const noexcept
  { return at_(key); }
  CelloView<T, 3> operator[] (std::size_t index) const noexcept
  { return at_(index); }

  /// Returns a reference to the mapped array associated with the specified
  /// index/key
  CelloView<T, 3> at(const std::string& key) const noexcept
  { return at_(key); }

  /// Checks whether the container holds the specified key
  bool contains(const std::string& key) const noexcept
  { return str_index_map_.contains(key); }

  /// Similar to `at`, but a slice of the array ommitting staled values is
  /// returned by value
  ///
  /// @note
  /// This is deprecated
  CelloView<T, 3> get(const std::string& key,
                      int stale_depth = 0) const noexcept
  { return get_(key, stale_depth); }

  /// Provided to help debug
  void print_summary() const noexcept;

  std::size_t size() const noexcept { return str_index_map_.size(); }

  /// Return the name of the instance (if it has one)
  const std::string& name() const noexcept {return name_;}

  /// Returns the length along a given dimension of each contained array
  ///
  /// @param dim Indicates the dimension for which we want the shape. This
  ///    should be 0, 1, or 2
  ///
  /// @note
  /// The program will fail if this method is invoked and the map holds zero
  /// elements.
  int array_shape(unsigned int dim) const noexcept
  { return arrays_.array_shape(dim); }

  /// Return a new map holding subsections of each array held by the map
  ///
  /// @param slc_z, slc_y, slc_x Instance of CSlice that specify that are
  ///    passed to the `subarray` method of each contained array.
  ViewMap<T> subarray_map(const CSlice &slc_z,
                          const CSlice &slc_y,
                          const CSlice &slc_x,
                          const std::string& name = "") const;

  /// Utility method offered for debugging purposes to check whether the
  /// key-order matches expectations.
  ///
  /// @param ref The list of reference keys
  /// @param raise_err When true, this will raise a fatal error if the order
  ///    doesn't match
  /// @param allow_smaller_ref When `true`, the key-order is only compared for
  ///    first `ref.size()`. When `false` (the default), the map must have the
  ///    same number of keys as `ref`.
  bool validate_key_order(const std::vector<std::string> &ref,
			  bool raise_err, bool allow_smaller_ref = false)
    const noexcept;

  /// Indicates if the contained arrays are stored in a single 4D array
  /// (Alternatively they can be stored as an array of pointers)
  bool contiguous_arrays() const noexcept { return arrays_.contiguous_items(); }

  /// Returns a shallow copy of the 4D array holding each contained array.
  ///
  /// The `n`th 3D subarray of `arrmap.get_backing_array()` is always a perfect
  /// alias of `arrmap[n]` (for any non-negative `n` less than `arrmap.size()`)
  ///
  /// @note
  /// The program will abort if this method is called on an object for which
  /// the `contiguous_arrays()` method returns `false`.
  CelloView<T, 4> get_backing_array() const noexcept
  { return arrays_.get_backing_array(); }

private: // helper methods

  /// This private constructor is used by subarray_map. It can skip some
  /// unnecessary work relating to initialization
  ViewMap(std::string name, const StringIndRdOnlyMap& str_index_map,
          ViewCollec<T>&& arrays)
    : name_(name),
      str_index_map_(str_index_map),
      arrays_(arrays)
  { validate_invariants_(); }

  void validate_invariants_() const noexcept;

  // The following 3 helper methods should NEVER be publically exposed
  CelloView<T, 3> at_(const std::string& key) const noexcept;
  CelloView<T, 3> at_(std::size_t index) const noexcept;
  CelloView<T, 3> get_(const std::string& key, int stale_depth)
    const noexcept;

private: // attributes
  // name_ is to help with debugging!
  std::string name_;

  // str_index_map_ maps the keys to the index
  StringIndRdOnlyMap str_index_map_;
  // arrays_ is the ordered collection of arrays_
  ViewCollec<T> arrays_;
};

//----------------------------------------------------------------------

template<typename T>
ViewMap<T>::ViewMap(std::string name, const std::vector<std::string> &keys,
                    const std::array<int,3>& shape)
  : name_(name),
    str_index_map_(keys),
    arrays_(keys.size(), shape)
{ }

//----------------------------------------------------------------------

template<typename T>
ViewMap<T>::ViewMap(std::string name, const std::vector<std::string> &keys,
                    const std::vector<CelloView<T,3>> &arrays)
  : name_(name),
    str_index_map_(keys),
    arrays_(arrays)
{
  ASSERT2("ViewMap::ViewMap",
          "keys and arrays have lengths %zu and %zu. They should be the same",
          (std::size_t)keys.size(), (std::size_t)arrays.size(),
          keys.size() == arrays.size());
}

//----------------------------------------------------------------------

template<typename T>
bool ViewMap<T>::validate_key_order(const std::vector<std::string> &ref,
                                    bool raise_err,
                                    bool allow_smaller_ref) const noexcept
{

  std::string name_string =
      (name_ == "") ? "" : (std::string(", \"") + name_ + std::string("\","));

  for (std::size_t i =0; i < std::min(ref.size(),str_index_map_.size()); i++){
    const std::string k = str_index_map_.key(i);
    bool equal = ref[i] == k;
    if (!equal && raise_err){
      print_summary();
      ERROR4("ViewMap::validate_key_order",
	     ("keys of ViewMap%s don't match expectations. At index %d, the "
	      "key is \"%s\" when it's expected to be \"%s\"\n"),
	     name_string.c_str(), (int)i, k.c_str(), ref[i].c_str());
    } else if (!equal){
      return false;
    }
  }

  if ((!allow_smaller_ref) && (ref.size() != str_index_map_.size())){
    if (raise_err){
      print_summary();
      ERROR3("ViewMap::validate_key_order",
	     "ViewMap%s doesn't have the expected number of keys. It has %d "
	     "keys. It's expected to have %d",
             name_string.c_str(), (int)ref.size(), (int)str_index_map_.size());
    } else {
      return false;
    }
  }
  return true;
}

//----------------------------------------------------------------------

template<typename T>
CelloView<T, 3> ViewMap<T>::at_(const std::string& key) const noexcept
{ return arrays_[str_index_map_.at(key)]; }

//----------------------------------------------------------------------

template<typename T>
CelloView<T, 3> ViewMap<T>::at_(const std::size_t index) const noexcept
{
  ASSERT("ViewMap::at_",
         "index must be less than or equal to the length of the ViewMap.",
         index < size());
  return arrays_[index];
}

//----------------------------------------------------------------------

template<typename T>
CelloView<T, 3> ViewMap<T>::get_(const std::string& key,
                                 int stale_depth) const noexcept
{
  ASSERT("ViewMap::get_", "stale_depth must be >= 0", stale_depth >= 0);
  CelloView<T, 3> view = this->at_(key);

  if (stale_depth > 0) {
    int min_len = stale_depth*2;
    ASSERT("ViewMap::get_","each dim of arr must exceed 2*stale_depth.",
           (min_len < view.shape(0)) && (min_len < view.shape(1)) &&
           (min_len < view.shape(2)));
    return view.subarray(CSlice(stale_depth, -stale_depth),
                         CSlice(stale_depth, -stale_depth),
                         CSlice(stale_depth, -stale_depth));
  } else {
    return view;
  }
}

//----------------------------------------------------------------------

template<typename T>
void ViewMap<T>::print_summary() const noexcept
{
  // TODO: would be nice to encode the template type in the string

  std::size_t my_size = size();
  if (my_size == 0){
    if (name_ == ""){
      CkPrintf("Nameless Empty View Map\n");
    } else {
      CkPrintf("\"%s\" Empty View Map\n", name_.c_str());
    }
    return;
  }

  if (name_ == ""){
    CkPrintf("Nameless View Map");
  } else {
    CkPrintf("\"%s\" View Map", name_.c_str());
  }

  CkPrintf(": entry_shape = (%d, %d, %d)\n{",
           array_shape(0), array_shape(1), array_shape(2));

  for ( std::size_t i = 0; i < my_size; i++){
    if (i != 0){
      CkPrintf(",\n ");
    }

    CelloView<T,3> array = arrays_[i];
    std::string key = str_index_map_.key(i);
    CkPrintf("\"%s\" : CelloView<T,3>(%p)", key.c_str(), (void*)array.data());
  }
  CkPrintf("}\n");
  fflush(stdout);
}

//----------------------------------------------------------------------

template<typename T>
ViewMap<T> ViewMap<T>::subarray_map(const CSlice &slc_z,
                                    const CSlice &slc_y,
                                    const CSlice &slc_x,
                                    const std::string& name) const
{
  return ViewMap<T>(name, str_index_map_,
                    arrays_.subarray_collec(slc_z, slc_y, slc_x));
}

//----------------------------------------------------------------------

template<typename T>
void ViewMap<T>::validate_invariants_() const noexcept
{
  // several invariants are alread enforced:
  // - StringIndRdOnlyMap implicitly enforces that there aren't any duplicate
  //   keys, and that a unique key is associated with each integer index from 0
  //   through (str_index_map_.size() - 1)
  // - ViewCollec enforces that all arrays have the same shape
  ASSERT("ViewMap::validate_invariants_",
         "str_index_map_ and arrays_ don't have the same length",
         arrays_.size() == str_index_map_.size());
}

#endif /* VIEW_VIEW_MAP_HPP */
