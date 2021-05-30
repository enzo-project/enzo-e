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

  const EFlt3DArray& operator[] (const std::string& key) const noexcept
  { return at(key); }

  const EFlt3DArray& operator[] (std::size_t index) const noexcept
  { return at(index); }

  const EFlt3DArray& at(const std::string& key) const noexcept;
  const EFlt3DArray& at(std::size_t index) const noexcept;

  bool contains(const std::string& key) const noexcept{
    return (str_index_map_.find(key) != str_index_map_.cend());
  }

  /// Similar to at, but a slice of the array ommitting staled values is
  /// returned by value
  EFlt3DArray get(const std::string& key,
                  int stale_depth = 0) const noexcept;

  /// Provided to help debug
  void print_summary() const noexcept;

  std::size_t size() const noexcept { return str_index_map_.size(); }

  /// Return the name of the instance (if it has one)
  const std::string& name() noexcept {return name_;}

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
