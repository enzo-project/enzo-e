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
};

//----------------------------------------------------------------------

/// @def      WRAPPED_LUT_INDEX_VALUE_T_SCALAR
/// @brief    Part of the WRAPPED_LUT_INDEX_VALUE group of macros that define
///           structs for forwarding for defining the enum members of
///           EnzoLUTWrapper for each of the actively advected quantities from
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

  enum quantites {
    #define ENTRY(name, math_type, category, if_advection)                 \
      WRAPPED_LUT_INDEX_VALUE_ ## if_advection ## _ ## math_type (name, LUT)
    FIELD_TABLE
    #undef ENTRY
  };

  // TODO: calculate the following automatically
  static constexpr std::size_t specific_start = LUT::specific_start;

  // TODO: calculate the following automatically
  static constexpr std::size_t NEQ=LUTIndexForward_::forward_NEQ<LUT>::value;

public: //associated static functions

  /// returns a set of integrable quantities included in the LUT
  static std::set<std::string> quantity_names(EnzoCenteredFieldRegistry &reg)
    noexcept;

  /// returns whether the LUT has any bfields
  static constexpr bool has_bfields(){
    return ( (quantites::bfield_i >= 0) ||
             (quantites::bfield_j >= 0) ||
             (quantites::bfield_k >= 0) );
  }
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

template<class LUT, class Function>
void unary_lut_for_each_(Function fn){
  #define ENTRY(name, math_type, category, if_advection)                      \
    LUT_UNARY_FUNC_ ## if_advection ## _ ## math_type (fn, LUT, name);
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
      if (index >= 0){
        if (reg.quantity_properties(name)){
          // current element is a SCALAR QUANTITY
          set.insert(name);
        } else if (reg.is_actively_advected_vector_component(name, true)){
          // current element is a VECTOR QUANTITY
          set.insert(name.substr(0,name.length()-2));
        }
      }
    };
  unary_lut_for_each_<EnzoRiemannLUTWrapper<LUT>>(fn);
  return set;
}

#endif /* ENZO_ENZO_RIEMANN_LUT_WRAPPER_HPP */
