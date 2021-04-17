// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoLazyNestedPassiveScalarFieldList.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sat January 9 2021
/// @brief    [\ref Enzo] Declaration of EnzoNestedPassiveScalarFieldList class

#include <pup_stl.h>

#ifndef ENZO_ENZO_NESTED_PASSIVE_SCALAR_FIELD_LIST_HPP
#define ENZO_ENZO_NESTED_PASSIVE_SCALAR_FIELD_LIST_HPP


class EnzoNestedPassiveScalarFieldList {
  /// @class    EnzoNestedPassiveScalarFieldList
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Functor that builds and manages the lifetime of
  ///           a nested passive scalar list
  ///
  /// The nested passive scalar list holds lists of passive scalar field names
  /// (or keys). The first sublist holds the names of all quantities that
  /// undergo normal passive advection. Subsequent lists hold sets of names for
  /// scalars whose values must sum to 1.0 (like species).
  ///
  /// Ideally, the list of lists would be directly built and managed by the
  /// Hydro integrator method object. However, the passive scalars will not all
  /// be known when the Hydro Integrator is initialized (because other method
  /// objects may be initialized later). Maybe this functionallity should be
  /// built into the Cello layer?
  ///
  /// @par Expected Usage
  /// An instance of this object is expected to be stored as a non-const member
  /// of a `Method` class. The `get_list` method should NOT be invoked in the
  /// constructor. The first invocation should during or after the first
  /// invocation of the `Method` class's `timestep` or `compute` method.
  ///
  /// @par
  /// The const variant of `get_list` (which is always used in `timestep`)
  /// will always build a new instance of the nested list, unless one has been
  /// cached during a call to non-const version of `get_list`.
  ///
  /// @par
  /// The underlying assumption in this object is that the collection of
  /// passive scalars does NOT change after the first time that they are
  /// retrieved. NOTE: It might be worth writing a debugger mode to check that
  /// this invariant is enforced.

public:
  EnzoNestedPassiveScalarFieldList()
    : initialized_(false),
      nested_names_(nullptr)
  { }

  /// Retrieve the nested passive scalar field list
  std::shared_ptr<const std::vector<str_vec_t>> get_list()
    noexcept
  {
    if (!initialized_){
      nested_names_ = build_nested_list_();
      initialized_ = true;
    }
    return nested_names_;
  }

  /// Retrieve the nested passive scalar field list
  ///
  /// This method is required to support the timestep method of Method objects.
  std::shared_ptr<const std::vector<str_vec_t>> get_list()
    const noexcept
  { return (initialized_) ? nested_names_ : build_nested_list_(); }

  void pup(PUP::er &p);

private:

  /// Helper method where the nested list is actually constructed
  static std::shared_ptr<const std::vector<str_vec_t>>
    build_nested_list_() noexcept;

private:
  bool initialized_;
  std::shared_ptr<const std::vector<str_vec_t>> nested_names_;
};


#endif /* ENZO_ENZO_NESTED_PASSIVE_SCALAR_FIELD_LIST_HPP */
