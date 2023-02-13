// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoEOSVariant.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Feb 10 2023
/// @brief    [\ref Enzo] Declaration of class EnzoEOSVariant

#ifndef ENZO_ENZO_EOSVARIANT_HPP
#define ENZO_ENZO_EOSVARIANT_HPP

#include <memory> // std::unique_ptr
#include <typeinfo> // needed for use with typeid
#include <type_traits> // std::add_pointer, std::is_pointer, std::remove_cv,
                       // std::is_same

class EnzoEOSVariant{

  /// @class    EnzoEOSVariant
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Represents a type-safe union of the EOS class
  ///
  /// This class is modelled after the std::variant template class and the
  /// associated std::visit template functions that is introduced to the
  /// standard library in C++ 17. When we eventually upgrade the codebase to
  /// use C++17, we should consider replacing this clas with those constructs
  /// directly.
  ///
  /// There are several advantages to using a union-based approach over an
  /// inheritance-based approach

public:

  EnzoEOSVariant& operator=(const EnzoEOSVariant&) = delete;
  EnzoEOSVariant& operator=(EnzoEOSVariant&&) = delete;

  /// return whether the ``EnzoEOSVariant`` holds an instance of type ``T``
  template <class T>
  bool holds_eos_type() const { return (get_if<T>() != nullptr); }

  /// If the ``EnzoEOSVariant`` holds an instance of type ``T``, returns a
  /// ``const`` reference to that instance. Otherwise the program aborts with
  /// an error.
  ///
  /// use the get_if method if you're uncertain about the contained type
  template <class T>
  const T& get() const{
    const T* temp = get_if<T>();
    if (temp == nullptr){
      // we store typeid(T).name() in a variable to avoid lifetime issues
      std::string type_name = typeid(T).name();
      ERROR1("EnzoEOSStructIdeal::get",
             "The EnzoEOSStructIdeal object doesn't contain an instance of %s",
             type_name.c_str());
    }
    return *temp;
  }

  /// If the ``EnzoEOSVariant`` holds an instance of type ``T``, returns a
  /// pointer to that instance (the pointer is readonly). Otherwise, returns a
  /// nullptr.
  ///
  /// If the type is unrecognized (i.e. ``EnzoEOSVariant`` is unable to hold
  /// an instance of that type), the program aborts with an error.
  ///
  /// @note
  /// The pointer's lifetime is linked to the ``EnzoEOSVariant``'s lifetime
  template <class T>
  std::add_pointer<const T>::type get_if() const
  {
    static_assert(false == std::is_pointer<T>::value,
                  "EnzoEOSVariant::get_if shouldn't be passed a pointer type");

    using cleanT = std::remove_cv<T>::type; // strip const qualification
    if (std::is_same<cleanT, EnzoEOSStructIdeal>::value){
      return ideal_eos_.get();
    } else if (std::is_same<cleanT, EnzoEOSStructIsothermal>::value){
      return isothermal_eos_.get();
    } else {
      // we store typeid(cleanT).name() in a variable to avoid lifetime issues
      std::string type_name = typeid(cleanT).name();
      ERROR1("EnzoEOSVariant::get_if", "The class is unaware of the type: %s",
             type_name.c_str());
    }
  }

private:
  // list of all of the different types of Equation Of State objects that can
  // be held. We currently use unique_ptrs so that we don't have to write a
  // destructor. All but one of these will be a nullptr

  std::unique_ptr<EnzoEOSStructIdeal> eos_ideal_;
  std::unique_ptr<EnzoEOSStructIsothermal> eos_isothermal_;

};

#endif /* ENZO_ENZO_EOSVARIANT_HPP */
