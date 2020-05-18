// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannUtils.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 15 2020
/// @brief    [\ref Enzo] Implementation of EnzoRiemannLUTWrapper

#ifndef ENZO_ENZO_RIEMANN_LUT_WRAPPER_HPP
#define ENZO_ENZO_RIEMANN_LUT_WRAPPER_HPP

//----------------------------------------------------------------------

/// @typedef lutarray
/// @brief   Specialization of std::array to be used to hold enzo_floats
///          associated with lookup tables (used with Riemann solvers)
template<class LUT>
using lutarray = std::array<enzo_float, LUT::NEQ>;

//----------------------------------------------------------------------

// Yields a combined token
#define COMBINE(prefix, suffix) prefix##suffix
#define COMBINE3(first, second, third) first##second##third
// Yields a stringized version of the combined token
#define STRINGIZE_COMBINE(prefix,suffix) STR_C_1(COMBINE(prefix,suffix))
#define STR_C_1(s) STR_C_2(s)
#define STR_C_2(s) #s

//----------------------------------------------------------------------

/// @def      CREATE_ENUM_VALUE_FORWARDER
/// @brief    Given an argument, {name}, this macro defines a struct called
///           forward_{name} that retrieves the value of the enum member
///           called {name} from a given struct/class. If the member is not
///           present, indicates a value of minus 1. This implementation is
///           loosely inspired by https://stackoverflow.com/a/16000226

#define CREATE_ENUM_VALUE_FORWARDER(name)                                     \
  /* Default case: target template doesn't have a member named {name}. */     \
  template <typename T, typename = int>                                       \
  struct COMBINE(forward_, name)                                              \
  {static constexpr int value = -1; };                                        \
                                                                              \
  /* Specialized case: */                                                     \
  template <typename T>                                                       \
  struct COMBINE(forward_, name) <T, decltype((void) T::name, 0)>             \
  { static constexpr int value = T::name; };

//----------------------------------------------------------------------
/// @def      LUT_INDEX_FORWARDER_T_SCALAR
/// @brief    Part of the LUT_INDEX_FORWARDER group of macros that define
///           structs for forwarding the values of enum members of a target
///           struct/class named for actively advected quantities from
///           FIELD_TABLE
#define LUT_INDEX_FORWARDER_T_SCALAR(name) CREATE_ENUM_VALUE_FORWARDER(name);
#define LUT_INDEX_FORWARDER_T_VECTOR(name)                                    \
  CREATE_ENUM_VALUE_FORWARDER(COMBINE(name, _i));                             \
  CREATE_ENUM_VALUE_FORWARDER(COMBINE(name, _j));                             \
  CREATE_ENUM_VALUE_FORWARDER(COMBINE(name, _k));
#define LUT_INDEX_FORWARDER_F_SCALAR(name) /* ... */
#define LUT_INDEX_FORWARDER_F_VECTOR(name) /* ... */

namespace LUTIndexForward_ {
  #define ENTRY(name, math_type, category, if_advection)                      \
    LUT_INDEX_FORWARDER_##if_advection ## _ ## math_type (name);
  FIELD_TABLE
  #undef ENTRY

  // forward number of equations
  CREATE_ENUM_VALUE_FORWARDER(NEQ);

  // forward the first index holding specific quantities
  CREATE_ENUM_VALUE_FORWARDER(specific_start);
};

//----------------------------------------------------------------------

/// @def      WRAPPED_LUT_INDEX_VALUE_T_SCALAR
/// @brief    Part of the WRAPPED_LUT_INDEX_VALUE group of macros that define
///           is used for defining that values of enum members in
///           `EnzoRiemannLUTWrapper` by forwarding the enum member values from
///           a separate LUT for each of the actively advected quantities from
///           FIELD_TABLE
#define WRAPPED_LUT_INDEX_VALUE_T_SCALAR(name, LUT)                           \
  name = LUTIndexForward_::COMBINE(forward_, name)<LUT>::value,
#define WRAPPED_LUT_INDEX_VALUE_T_VECTOR(name, LUT)                           \
  COMBINE(name,_i) = LUTIndexForward_::COMBINE3(forward_,name,_i)<LUT>::value,\
  COMBINE(name,_j) = LUTIndexForward_::COMBINE3(forward_,name,_j)<LUT>::value,\
  COMBINE(name,_k) = LUTIndexForward_::COMBINE3(forward_,name,_k)<LUT>::value,
#define WRAPPED_LUT_INDEX_VALUE_F_SCALAR(name, LUT) /* ... */
#define WRAPPED_LUT_INDEX_VALUE_F_VECTOR(name, LUT) /* ... */

//----------------------------------------------------------------------

template <typename InputLUT>
struct EnzoRiemannLUT{
  /// @class    EnzoRiemannLUT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides a compile-time lookup table that maps the
  ///           components of a subset of the actively advected quantities from
  ///           FIELD_TABLE to unique indices
  ///
  /// This is a template class that provides the following features at compile
  /// time:
  ///    - a lookup table (LUT) that maps the names of components of a subset
  ///      of the actively advected quantities defined in FIELD_TABLE to
  ///      unique, contiguous indices.
  ///    - the number of quantity components included in the table
  ///    - a way to iterate over just the conserved quantities or specific
  ///      quantities values that are stored in an array using these mapping
  ///    - a way to query which of the actively advected quantities in
  ///      FIELD_TABLE are not included in the LUT
  ///
  /// These feature are provided via the definition of publicly accessible
  /// integer constants in every specialization of the template class. All
  /// specializations have:
  ///    - a constant called `NEQ` equal to the number of quantity
  ///      components included in the lookup table
  ///    - a constant called `specific_start` equal to the number of components
  ///      of conserved quantities included in the lookup table
  ///    - `qkey` constants, which include constants named for the components
  ///      of ALL actively advected quantities in FIELD_TABLE. A constant
  ///      associated with a SCALAR quantity, `{qname}`, is simply called
  ///      `{qname}` while constants associated with a vector quantity
  ///      `{qname}` are called `{qname}_i`, `{qname}_j`, and `{qname}_k`.
  ///
  /// The `qkey` constants serve as both the keys of the lookup table and a
  /// way to check whether a component of an actively advected quantity is
  /// included in the table. Their values are satisfy the following conditions:
  ///    - All constants named for values corresponding to quantities included
  ///      in the table have values of `-1`
  ///    - All constants named for conserved quantities have unique integer
  ///      values in the internal `[0,specific_start)`
  ///    - All constants named for specific quantities have unique integer
  ///      values in the interval `[specific_start, NEQ)`
  ///
  /// The lookup table is always expected to include density and the 3 velocity
  /// components. Although it may not be strictly enforced (yet), the lookup
  /// table is also expected to include either all 3 components of a vector
  /// quantity or None of them.
  ///
  /// This template class also provides a handful of helpful static methods to
  /// programmatically probe the table's contents at runtime and validate that
  /// the above requirements are specified.
  ///
  /// @tparam InputLUT Type that defines the `NEQ` constant, the
  ///     `specific_start` constants and the `qkey` constants that correspond
  ///     to the actively advected quantities that are actually included in the
  ///     lookup table. All of the constant values will be directly reused by
  ///     the resulting template specialization and therefore must meet the
  ///     criteria defined above. Any undefined `qkey` constants are assumed to
  ///     not be included in the lookup table and will be defined within the
  ///     template specialization to have values of `-1`. This type can be
  ///     implemented with either an unscoped enum OR a class that includes
  ///     some combination of publicly accessible unscoped enums and static
  ///     constant member variables.
  ///
  /// @par Examples
  /// For the sake of the example, let's assume we have a type `MyInputLUT`
  /// that's defined as:
  /// @code
  ///      struct MyIntLUT {
  ///        enum vals { density=0, velocity_i, velocity_j, velocity_k,
  ///                    total_energy, NEQ, specific_start = 1};
  ///      };
  /// @endcode
  /// To access the index associated with density or the jth component of
  /// velocity, one would evaluate:
  /// @code
  ///      int density_index = EnzoRiemannLUT<MyInLUT>::density; //=0
  ///      int vj_index = EnzoRiemannLUT<MyInLUT>::velocity_j;   //=2
  /// @endcode
  ///
  /// @par
  /// It makes more sense to talk about the use of this template class when we
  /// have a companion array. For convenience, the alias template
  /// `lutarray<LUT>` type is defined. The type,
  /// `lutarray<EnzoRiemannLUT<InputLUT>>` is an alias of the type
  /// `std::array<enzo_float, EnzoRiemannLUT<InputLUT>::NEQ>;`.
  ///
  /// @par
  /// As an example, imagine that the total kinetic energy density needs to be
  /// computed at a single location from an values stored in an array, `prim`,
  /// of type `lutarray<EnzoRiemannLUT<MyInLUT>>`. The resulting code to
  /// do that might look something like:
  /// @code
  ///      using LUT = EnzoRiemannLUT<MyInLUT>;
  ///      enzo_float v2 = (prim[LUT::velocity_i] * prim[LUT::velocity_i] +
  ///                       prim[LUT::velocity_j] * prim[LUT::velocity_j] +
  ///                       prim[LUT::velocity_k] * prim[LUT::velocity_k]);
  ///      enzo_float kinetic = 0.5 * prim[LUT::density] * v2;
  /// @endcode
  ///
  /// @par
  /// This template class makes it very easy to write generic code that can be
  /// reused for multiple different lookup table by making the lookup table
  /// itself into a template argument. Consider the case where a single
  /// template function is desired to compute the total non-thermal energy
  /// density at a single location for an arbitrary lookup table. To do this,
  /// one might write the following code:
  /// @code
  ///      template <class LUT>
  ///      enzo_float calc_nonthermal_edens(lutarray<LUT> prim)
  ///      {
  ///        enzo_float v2 = (prim[LUT::velocity_i] * prim[LUT::velocity_i] +
  ///                         prim[LUT::velocity_j] * prim[LUT::velocity_j] +
  ///                         prim[LUT::velocity_k] * prim[LUT::velocity_k]);
  ///
  ///        enzo_float bi = (LUT::bfield_i >= 0) ? prim[LUT::bfield_i] : 0;
  ///        enzo_float bj = (LUT::bfield_j >= 0) ? prim[LUT::bfield_j] : 0;
  ///        enzo_float bk = (LUT::bfield_k >= 0) ? prim[LUT::bfield_k] : 0;
  ///        enzo_float b2 = bi*bi + bj*bj + bk*bk;
  ///
  ///        return 0.5(v2*prim[LUT::density] + b2);
  ///      }
  /// @endcode
  ///
  /// @par
  /// This final example will show the value of grouping the conserved and
  /// specific quantities. To implement some Riemann solvers, it is useful to
  /// have a generic, reusable function that constructs an array of all
  /// quantities that are in the lookup table in their conserved form. The
  /// following function was used to accomplish this at one point
  /// @code
  ///      template <class LUT>
  ///      lutarray<LUT> compute_conserved(const lutarray<LUT> prim) noexcept
  ///      {
  ///        lutarray<LUT> cons;
  ///        for (std::size_t i = 0; i < LUT::specific_start; i++) {
  ///          cons[i] = prim[i];
  ///        }
  ///        for (std::size_t i = LUT::specific_start; i < LUT::NEQ; i++) {
  ///          cons[i] = prim[i] * prim[LUT::density];
  ///        }
  ///        return cons;
  ///      }
  /// @endcode

public:
  /// initialize the lookup table entries for ever actively advected scalar
  /// quantity and component of a vector quantity in FIELD_TABLE by forwarding
  /// values from InputLUT (entries not included in InputLUT default to -1)
  enum qkey {
    #define ENTRY(name, math_type, category, if_advection)                 \
      WRAPPED_LUT_INDEX_VALUE_##if_advection##_##math_type (name, InputLUT)
    FIELD_TABLE
    #undef ENTRY
  };

  /// The entry with the minimum (non-negative) index corresponding to a
  /// specific quantity. All conserved quantity entries must have smaller
  /// indices (defaults to -1 if not explicitly specified in InputLUT)
  static constexpr std::size_t specific_start =
    LUTIndexForward_::forward_specific_start<InputLUT>::value;

  /// The total number of entries in the InputLUT (with non-negative indices).
  /// (defaults to -1 if not explicitly set in InputLUT)
  static constexpr std::size_t NEQ =
    LUTIndexForward_::forward_NEQ<InputLUT>::value;

  // perform some sanity checks:
  static_assert(qkey::density >= 0,
                "InputLUT must have an entry corresponding to density");
  static_assert((qkey::velocity_i >= 0 &&
                 qkey::velocity_j >= 0 &&
                 qkey::velocity_k >= 0),
                "InputLUT must have entries for each velocity component.");
  static_assert(specific_start > 0,
                "InputLUT::specific_start was not set to a positive value");
  static_assert(specific_start < NEQ,
                "InputLUT::NEQ was not set to a value exceeding "
                "InputLUT::specific_start");

private:
  /// This is not meant to be constructed.
  ///
  /// if it becomes necessary in the future to construct it in the future,
  /// there is no harm in doing so
  EnzoRiemannLUT() {}

public: //associated static functions

  /// returns a set of integrable quantities included in the InputLUT
  static std::set<std::string> quantity_names(EnzoCenteredFieldRegistry &reg)
    noexcept;

  /// returns whether the LUT has any bfields
  static constexpr bool has_bfields(){
    return (qkey::bfield_i>=0) || (qkey::bfield_j>=0) || (qkey::bfield_k>=0);
  }

  /// for each quantity in FIELD_TABLE, this passes the associated passes the
  /// associated name and index from the wrapped InputLUT to the function `fn`
  ///
  /// The function should expect (std::string name, int index). In cases where
  /// quantites in FIELD_TABLE do not appear in the wrapped InputLUT, an index
  /// of -1 is passed.
  template<class Function>
  static void for_each_entry(Function fn) noexcept;

  /// a function that performs a check to make sure that the InputLUT satisfies
  /// all assumptions. If it doesn't, an error is raised
  static void validate(EnzoCenteredFieldRegistry &reg) noexcept;
};

//----------------------------------------------------------------------

/// @def      LUT_UNARY_FUNC_T_SCALAR
/// @brief    Part of the LUT_UNARY_FUNC group of macros that apply a
///           unary function to the entries of named for advection related
///           quantities in FIELD_TABLE
#define LUT_UNARY_FUNC_T_SCALAR(func, LUT, name) func(#name, LUT::name)
#define LUT_UNARY_FUNC_T_VECTOR(func, LUT, name)                              \
  func(STRINGIZE_COMBINE(name,_i), LUT::COMBINE(name,_i));                    \
  func(STRINGIZE_COMBINE(name,_j), LUT::COMBINE(name,_j));                    \
  func(STRINGIZE_COMBINE(name,_k), LUT::COMBINE(name,_k))
#define LUT_UNARY_FUNC_F_SCALAR(func, LUT, name) /* ... */
#define LUT_UNARY_FUNC_F_VECTOR(func, LUT, name) /* ... */

//----------------------------------------------------------------------
/// template helper function that applies unary function on the enum members of
/// a lookup table that have been named for advection related quantities in
/// FIELD_TABLE
/// @param fn Unary function or that accepts the name of the  member and the
/// value of the member as arguments

template<class LUT>
template<class Function>
void EnzoRiemannLUT<LUT>::for_each_entry(Function fn) noexcept{
  #define ENTRY(name, math_type, category, if_advection)                      \
    LUT_UNARY_FUNC_##if_advection##_##math_type (fn,                          \
                                                 EnzoRiemannLUT<LUT>,  \
                                                 name);
  FIELD_TABLE
  #undef ENTRY
}

//----------------------------------------------------------------------

template <class InputLUT>
std::set<std::string> EnzoRiemannLUT<InputLUT>::quantity_names
(EnzoCenteredFieldRegistry &reg) noexcept
{
  std::set<std::string> set;

  auto fn = [&](std::string name, int index)
    {   
      if (index >= 0) {
        set.insert(reg.get_actively_advected_quantity_name(name, true));
      }
    };
  EnzoRiemannLUT<InputLUT>::for_each_entry(fn);
  return set;
}

//----------------------------------------------------------------------

template <class InputLUT>
void EnzoRiemannLUT<InputLUT>::validate(EnzoCenteredFieldRegistry &reg)
  noexcept
{ 
  // the elements in the array are default-initialized (they are each "")
  std::array<std::string, EnzoRiemannLUT<InputLUT>::NEQ> entry_names;

  // define a lambda function to execute for every member of lut
  auto fn = [&](std::string name, int index)
    {
      if ((index >= 0) && (index >= EnzoRiemannLUT<InputLUT>::NEQ)) {
        ERROR3("EnzoRiemannLUT<InputLUT>::validate",
               ("The value of %s, %d, is greater than or equal to "
                "InputLUT::NEQ, %d"),
               name.c_str(), index, EnzoRiemannLUT<InputLUT>::NEQ);
      } else if (index >= 0) {
        if (entry_names[index] != ""){
          ERROR3("EnzoRiemannLUT<InputLUT>::validate",
                 "%s and %s both have values of %d",
                 name.c_str(), entry_names[index], index);
        }
        entry_names[index] = name;
      }
    };
  

  EnzoRiemannLUT<InputLUT>::for_each_entry(fn);

  std::size_t max_conserved =  0;
  std::size_t min_specific =  EnzoRiemannLUT<InputLUT>::NEQ;

  for (std::size_t i = 0; i < entry_names.size(); i++){
    std::string name = entry_names[i];
    if (name == ""){
      ERROR2("EnzoRiemannLUT<LUT>::validate",
             "The value of NEQ, %d, is wrong. There is no entry for index %d",
             (int)EnzoRiemannLUT<InputLUT>::NEQ, (int)i);
    }
    // check that name != ''
    std::string quantity = reg.get_actively_advected_quantity_name(name, true);
    FieldCat category;
    reg.quantity_properties(quantity, NULL, &category, NULL);

    if ((i == 0) && (category != FieldCat::conserved)) {
      ERROR("EnzoRiemannLUT<InputLUT>::validate",
            ("the lookup table's entry for index 0 should correspond to a "
             "conserved quantity"));
    } else if (((i+1) == entry_names.size()) &&
               (category != FieldCat::specific)) {
      ERROR("EnzoRiemannLUT<InputLUT>::validate",
            ("the lookup table's entry for the index InputLUT::NEQ-1 should "
             "correspond to a specific quantity"));
    }

    switch(category){
      case FieldCat::conserved : { max_conserved = std::max(max_conserved, i);
                                   break;
      }
      case FieldCat::specific  : { min_specific  = std::min(min_specific,  i);
                                   break;
      }
      case FieldCat::other : {
        ERROR1("EnzoRiemannLUT<InputLUT>::validate",
               ("%s corresponds to a quantity with FieldCat::other. Quantities "
                "of this category should not be included in a lookup table."),
               name.c_str());
      }
    }
  }
   

  // The assumption is that all LUTs contain at least 1 conserved quantity
  // (nominally density) and at least 1 specific quantity (nominally velocity)
  if ((max_conserved+1) != min_specific){
    ERROR("EnzoRiemannLUT<InputLUT>::validate",
          ("InputLUT's entries corresponding to conserved quantities are not "
           "all grouped together at indices smaller than those corresponding "
           "to specific quantities"));
  } else if (min_specific != EnzoRiemannLUT<InputLUT>::specific_start) {
    ERROR2("EnzoRiemannLUT<InputLUT>::validate",
           "InputLUT's specfic_start value should be set to %d, not %d.",
           (int)min_specific, (int)EnzoRiemannLUT<InputLUT>::specific_start);
  }
}

#endif /* ENZO_ENZO_RIEMANN_LUT_WRAPPER_HPP */
