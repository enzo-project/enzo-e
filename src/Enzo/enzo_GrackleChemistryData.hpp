// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_GrackleChemistryData.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com) 
/// @date     Feb 6 2023
/// @brief    [\ref Enzo] Declaration of the GrackleChemistryData wrapper class
///
/// IMPORTANT: we are intentionally assuming that "grackle.h" is NOT visible to
/// this file and we are intentionally avoiding to conditionally alter the
/// declaration based on whether Grackle is being linked in-use.
///
/// By limiting checks about whether Grackle is being linked to the source
/// files, we will eventually be able to rebuild Enzo-E with/without Grackle
/// without recompiling the entire codebase

#ifndef ENZO_GRACKLE_CHEMISTRY_DATA_HPP
#define ENZO_GRACKLE_CHEMISTRY_DATA_HPP

// unclear whether extern "C" is needed here
extern "C" {
  typedef chemistry_data;
}

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
  ///   2. the str_field_map_ attribute manages the lifetime of all strings
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

  /// copy assignment operator (performs a deepcopy)
  GrackleChemistryData& operator= (const GrackleChemistryData&);

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// factory method that builds GrackleChemistryData from a Parameters object
  ///
  /// @param[in] p the titular Parameters object
  /// @param[in] ext_int_params, ext_dbl_params, ext_str_params Maps holding
  ///     parameter name value pairs that were taken from other places (e.g.
  ///     the adiabatic index
  /// @param[in] unrelated_params names of parameters that (may) occur within
  ///     the "Method:grackle:*" group that should be ignored!
  ///
  /// @note
  /// The abundance of parameters that are used to configure chemistry_data
  /// make it EXTREMELY easy to make a small mistake. For that reason, this
  /// function is VERY aggressive about reporting any unexpected parameters as
  /// an error.
  GrackleChemistryData static from_parameters
    (Parameters& p,
     std::unordered_map<std::string, int> ext_int_params,
     std::unordered_map<std::string, double> ext_dbl_params,
     std::unordered_map<std::string, std::string> ext_str_params,
     std::unordered_set<std::string> unrelated_params) noexcept;

public: // public interface (implementation is independent of Grackle details)

  /// default destructor
  ~GrackleChemistryData() = default;

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
  { ptr_.swap(other.ptr_); str_field_map_.swap(str_allocs_); }

  /// returns a pointer to the managed chemistry_data struct
  ///
  /// This is primarily intended to be used when calling a grackle function.
  /// The pointer is only valid during the GrackleChemistryData's lifetime
  inline chemistry_data* get_ptr() { return ptr_.get(); }
  inline const chemistry_data* get_ptr() const { return ptr_.get(); }

  /// query a field-value of the chemistry_data struct
  ///
  /// @note
  /// specializations of this method follow the class declaration
  template<class T>
  T get(const std::string& field) const {
    ERROR("GrackleChemistryData::get",
          "template parameter must be int, double, or std::string")
  }

  /// updates a field-values of the chemistry_data struct
  ///
  /// @note
  /// specializations of this method follow the class declaration
  template<class T>
  void set(const std::string& field, const T& value) {
    ERROR("GrackleChemistryData::set",
          "template parameter must be int, double, or std::string")
  }

private: // methods

  int get_int_(const std::string& field) const noexcept;
  double get_double_(const std::string& field) const noexcept;
  std::string get_string_(const std::string& field) const noexcept;

  void set_int_(const std::string& field, const int& value) noexcept;
  void set_double_(const std::string& field, const double& value) noexcept;
  void set_string_(const std::string& field, const std::string& value) noexcept;

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
template<> inline int GrackleChemistryData::get<>
  (const std::string& field) const noexcept { return get_int_(field); }
template<> inline double GrackleChemistryData::get<>
  (const std::string& field) const noexcept { return get_dbl_(field); }
template<> inline std::string GrackleChemistryData::get<>
  (const std::string& field) const noexcept { return get_str_(field); }

// specializations for set
template<> inline void GrackleChemistryData::set<>
  (const std::string& f, const int& v) noexcept {return set_int_(f, v);}
template<> inline void GrackleChemistryData::set<>
  (const std::string& f, const double& v) noexcept {return set_dbl_(f, v);}
template<> inline void GrackleChemistryData::set<>
  (const std::string& f, const std::string& v) noexcept {return set_str_(f, v);}


#endif /* ENZO_GRACKLE_CHEMISTRY_DATA_HPP */
