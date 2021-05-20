// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEFltArrayMap.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon Oct 26 2020
/// @brief    [\ref Enzo] Implementation of Enzo's EnzoEFltArrayMap which is
///           used to hold collections of EFlt3DArrays
///
/// This largely exists to help simplify the transition of EnzoMethodMHDVlct's
/// implementation from using Groupings to using Maps. It may not be optimal to
/// define this as it's own class
///
/// I'm starting to think that a more optimal approach may be to perform all
/// insertions in the constructor and to thereafter freeze the contained
/// arrays (the entries of the arrays could still be modified).
///   - The class could then provide the interface of a const std::map.
///   - It could also potentially enforce that all contained arrays have a
///     fixed shape
///
/// This choice, also facillitates further optimizations based on the capacity
/// to guarantee the ordering of arrays held within an instance. If the
/// ordering is based off the requirements of the Riemann Solver and the class
/// offers a method to return a (const?) vector of arrays with that ordering,
/// optimizations could be made to the Riemann Solver (as well as
/// EnzoReconstructor and EnzoIntegrableUpdater). This could be most
/// efficiently achieved by internally storing the contained array within a
/// vector and also storing some kind of mapping relating the keys to indices
/// the vector.
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

  EnzoEFltArrayMap(std::string name = "")
    : name_(name)
  { }

  EFlt3DArray& operator[] (const std::string& key){ return map_[key]; }
  EFlt3DArray& operator[] (std::string&& key){ return map_[key]; }

  const EFlt3DArray& at(const std::string& key) const noexcept;

  bool contains(const std::string& key) const noexcept{
    return (map_.find(key) != map_.cend());
  }

  /// Similar to at, but a slice of the array ommitting staled values is
  /// returned by value
  EFlt3DArray get(const std::string& key,
                  int stale_depth = 0) const noexcept;

  /// Provided to help debug
  void print_summary() const noexcept;

  std::size_t size() const noexcept { return map_.size(); }

private: // attributes
  // name_ is to help with debugging!
  std::string name_;
  std::map<std::string, EFlt3DArray> map_;
};

#endif /* ENZO_ENZO_EFLT_ARRAY_MAP_HPP */
