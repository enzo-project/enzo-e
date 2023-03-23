// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_GrackleChemistryData.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com) 
/// @date     Feb 6 2023
/// @brief    [\ref Enzo] Declaration of the GrackleChemistryData wrapper class
///
/// IMPORTANT: we are intentionally assuming that "grackle.h" is NOT visible to
/// this file and we are intentionally avoiding the practice of conditionally
/// altering the declaration based on whether Grackle is being linked.
///
/// By limiting checks about whether Grackle is being linked to the source
/// files, we will eventually be able to rebuild Enzo-E with/without Grackle
/// without recompiling the entire codebase

#ifndef ENZO_GRACKLE_CHEMISTRY_DATA_HPP
#define ENZO_GRACKLE_CHEMISTRY_DATA_HPP

// in the future, when the grackle header isn't included in the global header,
// we should always forward declare the chemistry data struct
#ifndef CONFIG_USE_GRACKLE
// unclear how necessary `extern "C"` is here
extern "C" { struct chemistry_data; };
#endif

class GrackleChemistryData {

  /// @class    GrackleChemistryData
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Wraps Grackle's chemistry_data struct. This provides a
  /// nice interface for accessing/modifying the struct's fields.
  ///
  /// Instances of this class satisfy 2 invariants:
  ///   1. the stored pointer to the ``chemistry_data`` struct is always valid;
  ///      it's not allowed to be a ``nullptr``. While lazily initializing the
  ///      pointer would make the move constructor cheaper, it's probably not a
  ///      worthwhile tradeoff (it would be REALLY easy to forget to lazily
  ///      initialize the pointer and then cause a bug)
  ///   2. the str_allocs_ attribute manages the lifetime of all strings
  ///      stored in fields of the ``chemistry_data`` struct. The only
  ///      exception is the default value (which may be a string literal or a
  ///      ``nullptr``)
  ///
  /// @note
  /// The implementation could definitely be optimized. In reality, the methods
  /// that are slow probably don't need to be much faster since they are mostly
  /// used during startup

public: // public interface (implementation is tied to Grackle details)

  /// construct a new GrackleChemistryData
  GrackleChemistryData();

  /// destroy a GrackleChemistryData instance
  ///
  /// this must be implemented in the source file where the full definition of
  /// chemistry_data is known (i.e. namely the fallback definition when Grackle
  /// isn't in use)
  ~GrackleChemistryData();

  /// copy assignment operator (performs a deepcopy)
  GrackleChemistryData& operator= (const GrackleChemistryData&);

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// factory method that builds GrackleChemistryData from a Parameters object
  ///
  /// @param[in] p the titular Parameters object
  /// @param[in] parameter_group string holding the full group name (usually it
  ///     will just be ``"Method:grackle"``.
  /// @param[in] forbid_leaf_names names of parameters that correspond to
  ///     grackle parameters that are not allowed to appear in the
  ///     "Method:grackle:*" group. These parameters are usually mutated after
  ///     the function is complete and takes the values from elsewhere.
  /// @param[in] ignore_leaf_names names of parameters that (may) occur within
  ///     the "Method:grackle:*" group that should be ignored!
  ///
  /// @note
  /// The abundance of parameters that are used to configure chemistry_data
  /// make it EXTREMELY easy to make a small mistake. For that reason, this
  /// function is VERY aggressive about reporting any unexpected parameters as
  /// an error.
  GrackleChemistryData static from_parameters
    (Parameters& p, const std::string& parameter_group,
     const std::unordered_set<std::string>& forbid_leaf_names,
     const std::unordered_set<std::string>& ignore_leaf_names) noexcept;

public: // public interface (implementation is independent of Grackle details)

  /// copy constructor (performs a deepcopy)
  ///
  /// @note
  /// this currently wraps the assignment operation
  GrackleChemistryData(const GrackleChemistryData& other)
    : GrackleChemistryData()
  { (*this) = other; }

  /// move constructor
  GrackleChemistryData(GrackleChemistryData&& other)
    : GrackleChemistryData()
  { this->swap(other); }

  /// assignment operation
  GrackleChemistryData& operator= (GrackleChemistryData&& other) noexcept
  { this->swap(other); return *this; }

  /// exchange the contents of this with other
  void swap(GrackleChemistryData& other) noexcept
  { ptr_.swap(other.ptr_); str_allocs_.swap(other.str_allocs_); }

  /// returns a pointer to the managed chemistry_data struct
  ///
  /// This is primarily intended to be used when calling a grackle function.
  /// The pointer is only valid during the GrackleChemistryData's lifetime
  inline chemistry_data* get_ptr() { return ptr_.get(); }
  inline const chemistry_data* get_ptr() const { return ptr_.get(); }

  /// query a field-value of the chemistry_data struct
  ///
  /// @note
  /// The choice to have this template method return by value and to define a
  /// separate template method for updating the values of parameters was made
  /// to provide a safe, uniform interface for managing parameters of all
  /// types. The alternative, having a single template method that returns a
  /// reference or pointer to the parameter, would generally produce issues for
  /// parameters of the string type.
  template<class T>
  T get(const std::string& field) const noexcept {
    std::pair<T,bool> tmp = this->try_get<T>(field);
    if (!tmp.second){
      std::string user_type = "";
      if (std::is_same<T, int>::value)         user_type = "int";
      if (std::is_same<T, double>::value)      user_type = "double";
      if (std::is_same<T, std::string>::value) user_type = "string";
      param_err_("GrackleChemistryData::get", user_type, field);
    }
    return tmp.first;
  }
  
  /// updates a field-values of the chemistry_data struct
  template<class T>
  void set(const std::string& field, const T& value) noexcept {
    if (!this->try_set<T>(field, value)){

      std::string user_type = "";
      if (std::is_same<T, int>::value)         user_type = "int";
      if (std::is_same<T, double>::value)      user_type = "double";
      if (std::is_same<T, std::string>::value) user_type = "string";

      param_err_("GrackleChemistryData::set", user_type, field);
    }
  }

  /// query a field-value of the chemistry_data struct. The first element holds
  /// a copy of the the value (upon success). The second element element is a
  /// boolean indicating whether the operation was successful
  ///
  /// @note
  /// The choice to have this template method return by value and to define a
  /// separate template method for updating the values of parameters was made
  /// to provide a safe, uniform interface for managing parameters of all
  /// types. The alternative, having a single template method that returns a
  /// reference or pointer to the parameter, would generally produce issues for
  /// parameters of the string type.
  ///
  /// @note
  /// specializations of this method follow the class declaration
  template<class T>
  std::pair<T,bool> try_get(const std::string& field) const noexcept {
    ERROR("GrackleChemistryData::try_get",
          "template parameter must be int, double, or std::string")
  }

  /// try to update a parameter value stored in the chemistry_data struct.
  ///
  /// @retval true the stored value was succesfully updated
  /// @retval false there is no know parameter of the specified type
  ///
  /// @note
  /// specializations of this method follow the class declaration
  template<class T>
  bool try_set(const std::string& field, const T& value) noexcept {
    ERROR("GrackleChemistryData::try_set",
          "template parameter must be int, double, or std::string")
  }

private: // methods

  /// static method used to generate nicely formatted error message when trying
  /// to access/modify a parameter that either doesn't have the user specified
  /// type or doesn't exist.
  ///
  /// Specifically, the error message reports whether the parameter
  /// exists at all or if the user-specified type is just wrong
  [[noreturn]] static void param_err_(const std::string& func_name,
                                      const std::string& user_type,
                                      const std::string& field);

  std::pair<int,bool> try_get_int_(const std::string& field) const noexcept;
  std::pair<double,bool> try_get_dbl_(const std::string& field) const noexcept;
  std::pair<std::string,bool> try_get_str_(const std::string& field)
    const noexcept;

  bool try_set_int_(const std::string& field, const int& value) noexcept;
  bool try_set_dbl_(const std::string& field, const double& value) noexcept;
  bool try_set_str_(const std::string& field,
                    const std::string& value) noexcept;

private: // attributes

  /// pointer to the chemistry_data struct
  std::unique_ptr<chemistry_data> ptr_;

  /// manages the lifetime of string fields used in chemistry_data
  ///
  /// We explicitly avoid using std::string because SSO (small string
  /// optimization) could cause all sorts of bugs in the future that would be
  /// VERY hard to debug
  ///
  /// The number of entries in this parameter is currently small. In the
  /// future, it may be sensible to replace vector with unordered_set
  std::vector<std::unique_ptr<char[]>> str_allocs_;
};

// specializations for get
template<> inline std::pair<int,bool> GrackleChemistryData::try_get<>
  (const std::string& field) const noexcept { return try_get_int_(field); }
template<> inline std::pair<double,bool> GrackleChemistryData::try_get<>
  (const std::string& field) const noexcept { return try_get_dbl_(field); }
template<> inline std::pair<std::string,bool> GrackleChemistryData::try_get<>
  (const std::string& field) const noexcept { return try_get_str_(field); }

// specializations for try_set
template<> inline bool GrackleChemistryData::try_set<>
  (const std::string& f, const int& v) noexcept
{return try_set_int_(f, v);}
template<> inline bool GrackleChemistryData::try_set<>
  (const std::string& f, const double& v) noexcept
{return try_set_dbl_(f, v);}
template<> inline bool GrackleChemistryData::try_set<>
  (const std::string& f, const std::string& v) noexcept
{return try_set_str_(f, v);}


#endif /* ENZO_GRACKLE_CHEMISTRY_DATA_HPP */
