// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEFltArrayMap.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon Oct 26 2020
/// @brief    [\ref Enzo] Implementation of Enzo's EnzoEFltArrayMap which is
///           used to hold collections of EFlt3DArrays
///
/// All insertions occur in the constructor, and the properties of the
/// contained arrays can't be mutated. This choice:
///   - facillitates enforcement that all contained arrays have a fixed shape
///   - makes it easier to order the entries in an arbitrary order. This could
///     lead to some optimizations in the Riemann Solver if the values are
///     initialized in the order expected by the Riemann Solver
///
/// If necessary, a number of optimizations could be made to the implementation
/// that might make key lookups faster. These optimizations could take
/// advantage of the following factors:
///    - Entries are never deleted and the contents won't be resized (if a hash
///      table can be resized, then the hash codes must stay the same or be
///      recomputed). These factors probably make a custom hash table using
///      open-addressing superior to std::map or std::unordered_map.
///    - Assumptions about the max key size and the max capacity of the map.
///      For example, if the max key size never exceeds ~22 characters and
///      there are never more than ~128 entries, it would probably be optimal
///      to store the strings in-place (improving cache locallity).
/// It would also be worth considering whether linear search is faster (since
/// the arrays are small.
///
/// To achieve similar results, instead of storing individual EFlt3DArrays
/// within vectors, one large instance of CelloArray<enzo_float,4> could be
/// stored. While that may help compiler optimizations, it may be too
/// restrictive.

#ifndef ENZO_ENZO_EFLT_ARRAY_MAP_HPP
#define ENZO_ENZO_EFLT_ARRAY_MAP_HPP

class EnzoEFltArrayMap {
  /// @class EnzoEFltArrayMap
  /// @ingroup Enzo
  /// @brief [\ref Enzo] Stores instances of EFlt3DArray

public: // interface

  EnzoEFltArrayMap(std::string name)
    : name_(name)
  { }

  EnzoEFltArrayMap()
    : EnzoEFltArrayMap("")
  { }

  /// Constructs a map that owns all of its own data.
  EnzoEFltArrayMap(std::string name, const std::vector<std::string> &keys,
                   const std::array<int,3>& shape);

  EnzoEFltArrayMap(const std::vector<std::string> &keys,
                   const std::array<int,3>& shape)
    : EnzoEFltArrayMap("", keys, shape)
  { }

  /// Constructs a map that wraps existing data.
  ///
  /// Each array must have the same shape.
  EnzoEFltArrayMap(std::string name, const std::vector<std::string> &keys,
                   const std::vector<EFlt3DArray> &arrays);

  EnzoEFltArrayMap(const std::vector<std::string> &keys,
                   const std::vector<EFlt3DArray> &arrays)
    : EnzoEFltArrayMap("", keys, arrays)
  { }

  /// Returns a reference to the mapped array associated with the specified
  /// index/key
  ///
  /// Unlike the analogous ``std::map::operator[]`` method, this cannot be used
  /// to insert a new value. This behaves similarly to ``std::map::at``
  CelloArray<const enzo_float, 3> operator[] (const std::string& key) const
    noexcept { return at_(key); }
  CelloArray<enzo_float, 3> operator[] (const std::string& key) noexcept
  { return at_(key); }
  CelloArray<const enzo_float, 3> operator[] (std::size_t index) const
    noexcept { return at_(index); }
  CelloArray<enzo_float, 3> operator[] (std::size_t index) noexcept
  { return at_(index); }

  /// Returns a reference to the mapped array associated with the specified
  /// index/key
  CelloArray<const enzo_float, 3> at(const std::string& key) const
    noexcept { return at_(key); }
  CelloArray<enzo_float, 3> at(const std::string& key) noexcept
    { return at_(key); }

  /// Checks whether the container holds the specified key
  bool contains(const std::string& key) const noexcept{
    return (str_index_map_.find(key) != str_index_map_.cend());
  }

  /// Similar to `at`, but a slice of the array ommitting staled values is
  /// returned by value
  CelloArray<enzo_float, 3> get(const std::string& key,
				int stale_depth = 0) noexcept
  { return get_(key, stale_depth); }
  CelloArray<const enzo_float, 3> get(const std::string& key,
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
  int array_shape(unsigned int dim) const noexcept;

  /// Return a new map holding subsections of each array held by the map
  ///
  /// @param slc_z, slc_y, slc_x Instance of CSlice that specify that are
  ///    passed to the `subarray` method of each contained array.
  EnzoEFltArrayMap subarray_map(const CSlice &slc_z,
                                const CSlice &slc_y,
                                const CSlice &slc_x,
                                const std::string& name = "");
  const EnzoEFltArrayMap subarray_map(const CSlice &slc_z,
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

private: // helper methods

  /// This private constructor is used by subarray_map. It can skip some
  /// unnecessary work relating to initialization
  EnzoEFltArrayMap(std::string name,
                   const std::map<std::string, unsigned int> &str_index_map,
                   const std::vector<std::string> &ordered_keys,
                   const std::vector<EFlt3DArray> &ordered_arrays)
    : name_(name),
      str_index_map_(str_index_map),
      keys_(ordered_keys),
      arrays_(ordered_arrays)
  { validate_invariants_(); }

  void validate_invariants_() const noexcept;

  // The following 3 helper methods should NEVER be publically exposed
  CelloArray<enzo_float, 3> at_(const std::string& key) const noexcept;
  CelloArray<enzo_float, 3> at_(std::size_t index) const noexcept;
  CelloArray<enzo_float, 3> get_(const std::string& key, int stale_depth)
    const noexcept;

private: // attributes
  // name_ is to help with debugging!
  std::string name_;
  
  // str_index_map_ maps the keys to the index
  std::map<std::string, unsigned int> str_index_map_;
  // keys_ is the ordered list of keys
  std::vector<std::string> keys_;
  // arrays_ is the ordered list of arrays_
  std::vector<EFlt3DArray> arrays_;
};

#endif /* ENZO_ENZO_EFLT_ARRAY_MAP_HPP */
