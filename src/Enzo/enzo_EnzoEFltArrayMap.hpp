// I'm not completely wedded to this idea of making a new class (as opposed to
// just using an existing class (like std::map)). But for now, this is the
// easiest approach to transition from Groupings to Maps

// Also not convinced that std::map is the best choice. I actually could
// imagine doing something like this with varying levels of complexity (e.g. we
// could wrap a vector of EFlt3DArray with some underlying order specified for
// the Riemann Solver, and use this class for source term calculations?)

// This does not have the standard library's container interface.

#ifndef ENZO_ENZO_EFLT_ARRAY_MAP_HPP
#define ENZO_ENZO_EFLT_ARRAY_MAP_HPP

class EnzoEFltArrayMap {
  /// @class EnzoEFltArrayMap
  /// @ingroup Enzo
  /// @brief [\ref Enzo] Stores instances of EFlt3DArray

public: // interface

  EFlt3DArray& operator[] (const std::string& key){ return map_[key]; }
  EFlt3DArray& operator[] (std::string&& key){ return map_[key]; }

  const EFlt3DArray& at(const std::string& key) const noexcept;

  bool contains(const std::string& key) const noexcept{
    return (map_.find(key) != map_.cend());
  }

  // This is a crutch for transitioning from groupings. Should get rid of this
  // eventually.
  const EFlt3DArray get(const std::string& key,
                        int stale_depth = 0) const noexcept;

  // This is also crutch for transitioning from groupings. Should get rid of
  // this eventually.
  static EnzoEFltArrayMap from_grouping(Block * block,
                                        Grouping& grouping,
                                        const std::vector<std::string>& groups,
                                        int dim = -1)
    noexcept;

private: // attributes
  std::map<std::string, EFlt3DArray> map_;
};

#endif /* ENZO_ENZO_EFLT_ARRAY_MAP_HPP */
