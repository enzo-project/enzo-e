// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEOSVariant.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-02-18
/// @brief    [\ref Enzo] Declaration of the EnzoEOSVariant class

#ifndef ENZO_ENZO_EOS_VARIANT_HPP
#define ENZO_ENZO_EOS_VARIANT_HPP

//----------------------------------------------------------------------

// forward declare EnzoEOSVariant
class EnzoEOSVariant;

//----------------------------------------------------------------------

/// function implementing the PUP operation
///
/// @note
/// This is not implemented as a method of EnzoEOSVariant to reflect the fact
/// that it won't be a method of std::variant either. It also ensures that we
/// don't rely upon implementation details of eos_variant (which ensures a
/// slightly more robust, if slightly slower) implementation.
///
/// @par
/// In the future, when we update this class's implementation to make use of
/// std::variant, we will probably need to rewrite the implementation in a more
/// manual way that doesn't use the macros. However, there's also a strong
/// possibility that we will have phased out pup routines by that time
void pup(PUP::er &p, EnzoEOSVariant& variant) noexcept;

//----------------------------------------------------------------------
  
// Next we define some tools to help implement EnzoEOSVariant::visit
namespace detail{

  /// This is an empty type that is used to keep track of a list of EOS Variant
  /// types (it is used to help implement EnzoEOSVariant::visit)
  template<class... Ts>
  struct EnzoEOSTypePack{ };

  //----------------------------------------------------------------------

  // The remaining stuff helps determines the type returned by the callable
  // that is passed to visit (i.e. the Visitor) AND it checks the assumption
  // that the return type is consistent independent of which EOS-Alternative
  // is actually passed.

  //----------------------------------------------------------------------

  /// Helper class that does most of the "heavy lifting". It determine the
  /// return type from passing each of the types ``Ts`` to a callable type,
  /// Visitor, separately.
  ///
  /// A compile-time error occurs if there is any variation in the return type/
  ///
  /// This ensures that Visitor (i) explicitly defines an overload of
  /// operator() for each ``T`` in ``Ts`` OR (ii) defines operator() as a
  /// function template that can accept each ``T`` in ``Ts``
  template<class Visitor, class... Ts>
  struct returnT_helper_; // not defined

  // recursive case
  template< class Visitor, class Head, class... Tail >
  struct returnT_helper_<Visitor, Head, Tail...> :
    returnT_helper_<Visitor, Tail...>
  {
    static_assert
    (std::is_same<typename returnT_helper_<Visitor,Head>::type,
                  typename returnT_helper_<Visitor, Tail...>::type>::value,
     "The visitor doesn't return the same type for all of the EOS types");
  };

  // base case
  template<class Visitor, class Tail>
  struct returnT_helper_<Visitor, Tail> {
    using type = typename std::result_of<Visitor(Tail)>::type;
  };

  //----------------------------------------------------------------------

  // down below, we're really just writing a nice interface for returnT_helper_

  //----------------------------------------------------------------------

  /// Determines the return type for passing any of the alternatives that the
  /// variant may contain to the Visitor.
  template<class Visitor, class Variant>
  struct visit_returnT; // not defined

  template<class Visitor, class... Ts>
  struct visit_returnT<Visitor, EnzoEOSTypePack<Ts...>>
    : returnT_helper_<Visitor, Ts...>
  { };
}

//----------------------------------------------------------------------

class EnzoEOSVariant {

  /// @class    EnzoEOSVariant
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Acts as a typesafe union for different that holds
  ///    an object encapsulating a caloric equation of state.
  ///
  /// This class has an interface that is modelled after std::variant.
  ///
  /// When we ultimately update the version of C++ that Enzo-E depends upon to
  /// C++17, we should use std::variant within the implementation of this class
  /// OR replace this class entirely with std::variant.

private:
  // lists of EOS alternatives (update these whenever you add a new alternative
  // type)

  // EOS_ALTERNATIVE_LIST is a list of the various alternative eos types, that
  // is used as part of the "Xmacro" technique. Each entry in the list:
  // - is passed a unique index for the type (this acts as the tag used
  //   internally to keep track of which alternative is currently stored by the
  //   variant
  // - is passed the type of the alternative.
#define EOS_ALTERNATIVE_LIST                 \
  EOS_ENTRY_(0, EnzoEOSIdeal)                \
  EOS_ENTRY_(1, EnzoEOSIsothermal)

  // this is a list used to help implement EnzoEOSVariant::visit
  using EOSTypes_ = detail::EnzoEOSTypePack<EnzoEOSIdeal, EnzoEOSIsothermal>;

private: // other types used internally

  union PtrUnion_{
    #define EOS_ENTRY_(alt_index, type) type* alt_ ## alt_index;
    EOS_ALTERNATIVE_LIST
    #undef EOS_ENTRY_
  };

  // this is just a forward declaration of a type that we define at the end
  // of this Header (that's used to implement get_if
  template< typename T> struct access_helper_;

private:

  void cleanup_(){
    switch(tag_) {
      #define EOS_ENTRY_(IND, _) case IND: {delete ptr_union_.alt_##IND; break;}
      EOS_ALTERNATIVE_LIST
      #undef EOS_ENTRY_

      default: ERROR1("EnzoEOSVariant::cleanup", "unhandled tag: %d", tag_);
    }
  }

public: // interface

  // define converting constructors
  //   EnzoEOSVariant(EnzoEOSIdeal eos);
  //   EnzoEOSVariant(EnzoEOSIsothermal eos);
  //   ...
#define EOS_ENTRY_(IND, TYPE)                                                 \
  EnzoEOSVariant(TYPE eos)                                                    \
  {ptr_union_.alt_##IND = new TYPE {eos}; tag_=IND;}

  EOS_ALTERNATIVE_LIST
#undef EOS_ENTRY_

  /// Default constructor.
  ///
  /// We're choosing to create a valid isothermal eos. In the future,
  /// we probably don't want to make guarantees the contained eos is
  /// necessarily configured with sensible values...
  EnzoEOSVariant() : EnzoEOSVariant(EnzoEOSIsothermal{}) { }

  /// Copy constructor
  ///
  /// This could definitely be better optimized
  EnzoEOSVariant(const EnzoEOSVariant& other)
    : EnzoEOSVariant()
  {
    cleanup_(); // if we skipped the default constructor, we could also
                // skip this, but let's just be safe
    switch(other.tag_) {
      #define EOS_ENTRY_(IND, T)                                             \
        case IND: {                                                          \
          this->ptr_union_.alt_##IND = new T{*other.ptr_union_.alt_##IND};   \
          break;                                                             \
        }
      EOS_ALTERNATIVE_LIST
      #undef EOS_ENTRY_

      default: ERROR1("EnzoEOSVariant::EnzoEOSVariant",
                      "unhandled tag: %d", tag_);
    }

    this->tag_ = other.tag_;
  }

  /// copy assignment
  EnzoEOSVariant& operator=(const EnzoEOSVariant& other) {
    EnzoEOSVariant copy_of_other(other);
    this->swap(copy_of_other);
    return *this;
  }

  /// move constructor
  EnzoEOSVariant(EnzoEOSVariant&& other)
    : EnzoEOSVariant()
  { this->swap(other); }

  /// move assignment
  EnzoEOSVariant& operator=(EnzoEOSVariant&& other) {
    this->swap(other);
    return *this;
  }

  /// destructor
  ~EnzoEOSVariant() noexcept { cleanup_(); }

  /// swap contents with other
  void swap(EnzoEOSVariant& other) {
    std::swap(this->ptr_union_, other.ptr_union_);
    std::swap(this->tag_, other.tag_);
  }

  /// returns a pointer to the contained EOS object if the EOS object's type
  /// matches the specified template parameter. Otherwise, a ``nullptr`` is
  /// returned
  ///
  /// The program aborts with an error, if a template parameter is passed that
  /// doesn't correspond to a known type the ``this`` can EVER hold.
  ///
  /// @note
  /// This effectively backports C++17's std::get_if template function.
  template <typename T>
  T* get_if() noexcept;

  template <typename T>
  const T* get_if() const noexcept
  { return const_cast<EnzoEOSVariant*>(this)->get_if<T>(); }

  /// accessor method that returns a reference to the contained EOS object, if
  /// ``this`` curently holds the EOS object of type ``T``. Otherwise, the
  /// program aborts with an error message.
  ///
  /// @tparam T the type of the EOS object that is being accessed
  ///
  /// @note
  /// This acts as a backports for one of C++17's ``std::get`` template
  /// function that is implemented for ``std::variant``.
  template<typename T>
  T& get() noexcept
  {
    T* out = this->get_if<T>();
    if (out == nullptr) {
      ERROR("EnzoEOSVariant::get",
            "The object doesn't hold a reference to the specified type.");
    }
    return *out;
  }

  template <typename T>
  const T& get() const noexcept
  { return const_cast<EnzoEOSVariant*>(this)->get<T>(); }

  /// Checks whether ``this`` currently holds the alternative EOS type, T.
  /// If T is not a type that can be held, the program aborts with an error.
  ///
  /// @note
  /// This acts as a backports for one of C++17's ``std::holds_alternative``
  template<typename T>
  bool holds_alternative() const noexcept { return get_if<T>() != nullptr; }

  /// invokes the callable visitor, ``vis`` by passing the EOS instance held by
  /// ``*this``
  ///
  /// The visitor must accept any of the EOS variants passed as an argument, by
  /// value, and return an output with a consistent type for all of them.
  ///
  /// This is a somewhat crude backport of ``std::visit`` from C++17.
  ///
  /// @note
  /// We could alternatively allow the visitor to accept a const reference to
  /// the contained eos object. However, we expressly forbid mutation of the
  /// constained eos object.
  template<class Visitor>
  typename detail::visit_returnT<Visitor, EOSTypes_>::type visit(Visitor&& vis)
    const noexcept
  {
    // the following switch statement expands into:
    // switch(tag_) {
    //   case 0: {
    //     const EnzoEOSIdeal& ref = ptr_union_.alt_0; return vis(ref);
    //   }
    //   case 1: {
    //     const EnzoEOSIsothermal& ref = ptr_union_.alt_1; return vis(ref);
    //   }
    //   ...
    //   default: ERROR1(...);
    // }

    switch(tag_) {
      #define EOS_ENTRY_(IND, T)                                             \
        case IND: { const T& ref = *(ptr_union_.alt_##IND); return vis(ref); }
      EOS_ALTERNATIVE_LIST
      #undef EOS_ENTRY_

      default: ERROR1("EnzoEOSVariant::visit", "unhandled tag: %d", tag_);
    }
  }

private:

  // the fact that this class currently holds pointers to the different variant
  // types is an implementation detail. We need to implement the class's
  // semantics as though we're holding a union of value types
  PtrUnion_ ptr_union_;

  int tag_;
};

// down below we implement the machinery needed for EnzoEOSVariant::get_if()
//
// Invalid Naive Solution
// ======================
// naively one might try to implement ``EnzoEOSVariant::get_if()`` with
// something like the following:
//
// template <typename T>
// T* EnzoEOSVariant::get_if() noexcept {
//   switch (tag_){
//     case 0: return (std::is_same<T, EnzoEOSIdeal>::value)
//       ? ptr_union_.alt_0 : nullptr;
//     case 1: return (std::is_same<T, EnzoEOSIsothermal>::value)
//       ? ptr_union_.alt_1 : nullptr;
//     ...
//     default: ERROR1("EnzoEOSVariant::get_if()", "invalid tag %d", tag_);
//   }
// }
//
// However, the compiler will complain about how mismatching of typing between
// members of the union and T*. For example, if T = EnzoEOSIdeal, then we'll
// the compiler will complain that one of the branches tries to return
// ptr_union_.alt_1, which can't be coerced to ``EnzoEOSIdeal*``.
//
// Actual Solution
// ===============
// Instead, we implement ``EnzoEOSVariant::get_if()`` using full template
// specialization of the ``EnzoEOSVariant::access_helper_`` class template.
// - This class template just defines a single static function: ``get_ptr``.
// - As an aside, we MUST define specializations of a class template instead
//   of a stand-alone function template called ``get_ptr``, because the only
//   difference between different specializations versions of ``get_ptr`` is
//   the return-type.


// This is the primary definition of the EnzoEOSVariant::access_helper_ class
// template. We will define separate specializations for all known EOS types
// that EnzoEOSVariant can hold. Thus, this primary definition is only ever
// used if the user tries to access an unrecognized eos type
template <typename T>
struct EnzoEOSVariant::access_helper_ {
  static T* get_ptr(const EnzoEOSVariant::PtrUnion_& ptr_union, int tag)
  {ERROR("EnzoEOSVariant::access_helper_::get_ptr", "unrecognized eos type");}
};

// Here, we use macros to define a separate specialization of the
// EnzoEOSVariant::access_helper_ class template for each known EOS types that
// EnzoEOSVariant can hold.
#define EOS_ENTRY_(IND, TYPE)                                                 \
  template <>                                                                 \
  struct EnzoEOSVariant::access_helper_<TYPE> {                               \
    static TYPE* get_ptr(const EnzoEOSVariant::PtrUnion_& ptr_union, int tag) \
    { return (tag == IND) ? ptr_union.alt_##IND : nullptr; }                  \
  };

EOS_ALTERNATIVE_LIST
#undef EOS_ENTRY_

// finally, we use EnzoEOSVariant::access_helper_ to implement get_if
template <typename T>
T* EnzoEOSVariant::get_if() noexcept
{ return EnzoEOSVariant::access_helper_<T>::get_ptr(ptr_union_, tag_); }

#endif /* ENZO_ENZO_EOS_VARIANT_HPP */
