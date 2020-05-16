// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannUtils.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 15 2020
/// @brief    [\ref Enzo] Implementation of EnzoRiemannLUTWrapper

#ifndef ENZO_ENZO_RIEMANN_LUT_WRAPPER_HPP
#define ENZO_ENZO_RIEMANN_LUT_WRAPPER_HPP

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

template <typename LUT>
struct EnzoRiemannLUTWrapper{
  /// @class    EnzoRiemannLUTWrapper
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Wraps LUT and provides a compile-time LUT that
  ///           maps indices to each actively advected quantity in
  ///           FIELD_TABLE
  ///
  /// This a bare-bones template class with enum members named for different
  /// advection quantites listed in FIELD_TABLE that each correspond to a
  /// unique index. The mathematical type of the quantity determines the name
  /// of the member. A SCALAR simply corresponds to a member with a name
  /// matching the value in column 1. A VECTOR corresponds members named
  /// {column_1}_i, {column_1}_j, {column_2}_k.
  ///
  /// This template class wraps existing LUTs that just 
  ///
  /// As an example, to check if a lookup table, LUT, contains the density
  /// quantity evaluate
  /// @code
  ///      EnzoLUTQuantityInfo<LUT>::density;
  /// @endcode
  /// Likewise, to check if LUT contains the velocity quantity, evaluate
  /// @code
  ///      EnzoLUTQuantityInfo<LUT>::velocity;
  /// @endcode

public:
  /// initialize the lookup table entries for ever actively advected scalar
  /// quantity and component of a vector quantity in FIELD_TABLE by forwarding
  /// values from LUT (entries not included in LUT default to -1)
  enum quantities {
    #define ENTRY(name, math_type, category, if_advection)                 \
      WRAPPED_LUT_INDEX_VALUE_ ## if_advection ## _ ## math_type (name, LUT)
    FIELD_TABLE
    #undef ENTRY
  };

  /// The entry with the minimum (non-negative) index corresponding to a
  /// specific quantity. All conserved quantity entries must have smaller
  /// indices (defaults to -1 if not explicitly specified in LUT)
  static constexpr std::size_t specific_start =
    LUTIndexForward_::forward_specific_start<LUT>::value;

  /// The total number of entries in the LUT (with non-negative indices).
  /// (defaults to -1 if not explicitly set in LUT)
  static constexpr std::size_t NEQ = LUTIndexForward_::forward_NEQ<LUT>::value;

  // perform some sanity checks:
  static_assert(quantities::density >= 0,
                "LUT must have an entry corresponding to density");
  static_assert((quantities::velocity_i >= 0 &&
                 quantities::velocity_j >= 0 &&
                 quantities::velocity_k >= 0),
                "LUT must have entries for each velocity component.");
  static_assert(specific_start > 0,
                "specific_start was not set to a be positive value in LUT");
  static_assert(specific_start < NEQ,
                "NEQ was not set to a value exceeding specific_start in LUT");

public: //associated static functions

  /// returns a set of integrable quantities included in the LUT
  static std::set<std::string> quantity_names(EnzoCenteredFieldRegistry &reg)
    noexcept;

  /// returns whether the LUT has any bfields
  static constexpr bool has_bfields(){
    return ( (quantities::bfield_i >= 0) ||
             (quantities::bfield_j >= 0) ||
             (quantities::bfield_k >= 0) );
  }

  /// for each quantity in FIELD_TABLE, this passes the associated passes the
  /// associated name and index from the wrapped LUT to the function `fn`
  ///
  /// The function should expect (std::string name, int index). In cases where
  /// quantites in FIELD_TABLE do not appear in the wrapped LUT, an index of -1
  /// is passed.
  template<class Function>
  static void for_each_entry(Function fn) noexcept;

  /// a function that performs a check to make sure that the LUT satisfies all
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
void EnzoRiemannLUTWrapper<LUT>::for_each_entry(Function fn) noexcept{
  #define ENTRY(name, math_type, category, if_advection)                      \
    LUT_UNARY_FUNC_##if_advection##_##math_type (fn,                          \
                                                 EnzoRiemannLUTWrapper<LUT>,  \
                                                 name);
  FIELD_TABLE
  #undef ENTRY
}

//----------------------------------------------------------------------

template <class LUT>
std::set<std::string> EnzoRiemannLUTWrapper<LUT>::quantity_names
(EnzoCenteredFieldRegistry &reg) noexcept
{
  std::set<std::string> set;

  auto fn = [&](std::string name, int index)
    {   
      if (index >= 0) {
        set.insert(reg.get_actively_advected_quantity_name(name, true));
      }
    };
  EnzoRiemannLUTWrapper<LUT>::for_each_entry(fn);
  return set;
}

//----------------------------------------------------------------------

template <class LUT>
void EnzoRiemannLUTWrapper<LUT>::validate(EnzoCenteredFieldRegistry &reg)
  noexcept
{ 
  // the elements in the array are default-initialized (they are each "")
  std::array<std::string, EnzoRiemannLUTWrapper<LUT>::NEQ> entry_names;

  // define a lambda function to execute for every member of lut
  auto fn = [&](std::string name, int index)
    {
      if ((index >= 0) && (index >= EnzoRiemannLUTWrapper<LUT>::NEQ)) {
        ERROR3("EnzoRiemannLUTWrapper<LUT>::validate",
               "The value of %s, %d, is greater than or equal to LUT::NEQ, %d",
               name.c_str(), index, EnzoRiemannLUTWrapper<LUT>::NEQ);
      } else if (index >= 0) {
        if (entry_names[index] != ""){
          ERROR3("EnzoRiemannLUTWrapper<LUT>::validate",
                 "%s and %s both have values of %d",
                 name.c_str(), entry_names[index], index);
        }
        entry_names[index] = name;
      }
    };
  

  EnzoRiemannLUTWrapper<LUT>::for_each_entry(fn);

  std::size_t max_conserved =  0;
  std::size_t min_specific =  EnzoRiemannLUTWrapper<LUT>::NEQ;

  for (std::size_t i = 0; i < entry_names.size(); i++){
    std::string name = entry_names[i];
    if (name == ""){
      ERROR2("EnzoRiemannLUTWrapper<LUT>::validate",
             "The value of NEQ, %d, is wrong. There is no entry for index %d",
             (int)EnzoRiemannLUTWrapper<LUT>::NEQ, (int)i);
    }
    // check that name != ''
    std::string quantity = reg.get_actively_advected_quantity_name(name, true);
    FieldCat category;
    reg.quantity_properties(quantity, NULL, &category, NULL);

    if ((i == 0) && (category != FieldCat::conserved)) {
      ERROR("EnzoRiemannLUTWrapper<LUT>::validate",
            ("the lookup table's entry for index 0 should correspond to a "
             "conserved quantity"));
    } else if (((i+1) == entry_names.size()) &&
               (category != FieldCat::specific)) {
      ERROR("EnzoRiemannLUTWrapper<LUT>::validate",
            ("the lookup table's entry for the index LUT::NEQ-1 should "
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
        ERROR1("EnzoRiemannLUTWrapper<LUT>::validate",
               ("%s corresponds to a quantity with FieldCat::other. Quantities "
                "of this category should not be included in a lookup table."),
               name.c_str());
      }
    }
  }
   

  // The assumption is that all LUTs contain at least 1 conserved quantity
  // (nominally density) and at least 1 specific quantity (nominally velocity)
  if ((max_conserved+1) != min_specific){
    ERROR("EnzoRiemannLUTWrapper<LUT>::validate",
          ("The LUT's entries corresponding to conserved quantities are not "
           "all grouped together at indices smaller than those corresponding "
           "to specific quantities"));
  } else if (min_specific != EnzoRiemannLUTWrapper<LUT>::specific_start) {
    ERROR2("EnzoRiemannLUTWrapper<LUT>::validate",
           "The LUT's specfic_start value should be set to %d, not %d.",
           (int)min_specific, (int)EnzoRiemannLUTWrapper<LUT>::specific_start);
  }
}

#endif /* ENZO_ENZO_RIEMANN_LUT_WRAPPER_HPP */
