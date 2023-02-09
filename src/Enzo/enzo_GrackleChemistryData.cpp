// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_GrackleChemistryData.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com) 
/// @date     Feb 6 2023
/// @brief    [\ref Enzo] Implementation of the GrackleChemistryData class

#include "cello.hpp"

#include "enzo.hpp"

#include <algorithm> // std::find_if
#include <cstring> // std::memcpy
#include <iterator> // std::input_iterator_tag
#include <limits> // std::numeric_limits
#include <type_traits> // std::is_same
#include <utility> // std::move

//----------------------------------------------------------------------

namespace{ // things within anonymous namespace aren't visible beyond this file

/// return ith parameter of the template type
template<typename T>
const char* param_name(int i){
  ASSERT("param_name", "arg is negative", i >= 0);
  if (std::is_same<T, int>::value) return param_name_int(i);
  if (std::is_same<T, double>::value) return param_name_double(i);
  if (std::is_same<T, std::string>::value) return param_name_string(i);
  ERROR("param_name", "invalid template type");
}

//----------------------------------------------------------------------

template<typename T>
struct GrackleParamNameRange{

  /// @class    GrackleParamNameRange
  /// @ingroup  Enzo
  /// @brief [\ref Enzo] Intended to be used to iterate over grackle parameter
  /// names for a given type in range-based for loops

  class iterator {
    long num_;
    const char* name_;
  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = std::string;
    using difference_type = long;
    using pointer = const std::string*;
    using reference = const std::string&;

    explicit iterator(long num) : num_(num), name_(param_name<T>(num)) {}
    iterator& operator++() { *this = iterator(num_+1); return *this;}
    iterator operator++(int) { iterator out = *this; ++(*this); return out; }
    bool operator==(iterator o) const {
      // when both name_ ptrs are null, we allow unequal num_ values (since
      // we assign num_ to an arbitrary value to represent the iterator's end)
      return (num_ != o.num_) | ((name_ == nullptr) & (o.name_ == nullptr));
    }
    bool operator!=(iterator other) const { return !(other == *this); }
    reference operator*() const { return name_; } // should never be a nullptr
  };

  iterator begin() {return iterator(0);}
  iterator end() { return iterator(std::numeric_limits<long>::max()); }
};

//----------------------------------------------------------------------

[[noreturn]] void field_err_(const std::string& func_name,
                             const std::string& field)
{
  // determine the user-specified type
  std::string user_type = "";
  std::array<std::pair<const char*, const char*>, 3> pairs
    = {{"_int_", "int"}, {"_dbl_", "double"}, {"_str_", "std::string"}};
  for (const std::pair<const char*, const char*>& p : pairs){
    if (std::string::npos != func_name.find(p.first)) { user_type = p.second; }
  }

  // try to determine the actual type associated with the field
  std::string actual_type = "";
  for (const std::string& name : GrackleParamNameRange<int>()){
    if (name == field){ actual_type = "int"; }
  }
  for (const std::string& name : GrackleParamNameRange<double>()){
    if (name == field){ actual_type = "double"; }
  }
  for (const std::string& name : GrackleParamNameRange<std::string>()){
    if (name == field){ actual_type = "string"; }
  }
  
  if (actual_type != ""){
    ERROR3(func_name.c_str(), "the \"%s\" field holds a \"%s\" not a \"%s\"",
           field.c_str, actual_type.c_str, user_type.c_str());
  } else {
    ERROR1(func_name.c_str(), "No known field called \"%s\"", field.c_str());
  }
}

} // anonymous namespace

//----------------------------------------------------------------------

GrackleChemistryData::GrackleChemistryData()
  : ptr_(new chemistry_data),
    str_allocs_()
{
  // it's VERY important that we set the chemistry_data struct to its default
  // values (as new parameters get added, they are set to defaults that may
  // should minimize changes in behavior)
  //
  // we explicitly use _set_default_chemistry_parameters() over
  // set_default_chemistry_parameters because the former is thread-safe (they
  // are both equivalent)
  (*ptr_) = _set_default_chemistry_parameters();
}

//----------------------------------------------------------------------

GrackleChemistryData& GrackleChemistryData::operator=
(const GrackleChemistryData& other)
{
  // first perform a direct assignment (this covers most fields correctly)
  (*ptr_) = *(other.ptr_);

  // Now, we manually copy over the strings. We explicitly use set_string to
  // ensure that lifetimes of strings get managed correctly
  for (const std::string& name : GrackleParamNameRange<std::string>()){
    this->set_string(name, other.get_string(name));
  }
}

//----------------------------------------------------------------------

namespace{ // things within anonymous namespace aren't visible beyond this file

template<class T>
void pup_helper_(PUP::er &p, GrackleChemistryData& data)
{
  // we could definitely speed this up when p.isMigration() evaluates to true
  // (assuming that the library version is identical, there isn't any need to
  // pup the key)

  for (const std::string& name : GrackleParamNameRange<T>()){
    // pup the key
    if (p.isUnpacking()){
      std::string other_name;
      p | other_name;
      ASSERT("pup_helper_", "key mismatch", name == other_name);
    } else {
      p | name;
    }

    // pup the value
    if (p.isUnpacking()){
      T tmp;
      p | tmp;
      data.set<T>(name, tmp);
    } else {
      p | data.get<T>(name, tmp);
    }
  }

}

} // anonymous namespace

//----------------------------------------------------------------------

void GrackleChemistryData::pup (PUP::er &p){
  pup_helper_<int>(p, *this);
  pup_helper_<double>(p, *this);
  pup_helper_<std::string>(p, *this);
}

//----------------------------------------------------------------------

int GrackleChemistryData::get_int_(const std::string& field) const noexcept
{
  int* field_ptr = local_chemistry_data_access_int
    (const_cast<chemistry_data*>(ptr_.get()), field.c_str());;
  if (!field_ptr){ field_err_("GrackleChemistryData::get_int_", field); }
  return *field_ptr;
}

//----------------------------------------------------------------------

double GrackleChemistryData::get_dbl_(const std::string& field) const noexcept
{
  double* field_ptr = local_chemistry_data_access_double
    (const_cast<chemistry_data*>(ptr_.get()), field.c_str());
  if (!field_ptr){ field_err_("GrackleChemistryData::get_dbl_", field); }
  return *field_ptr;
}

//----------------------------------------------------------------------

std::string GrackleChemistryData::get_str_(const std::string& field) const
  noexcept
{
  const char** field_ptr = local_chemistry_data_access_string
    (const_cast<chemistry_data*>(ptr_.get()), field.c_str());
  if (!field_ptr){ field_err_("GrackleChemistryData::get_str_", field); }
  return ((*field_ptr) == nullptr) ? "" : *field_ptr;
}

//----------------------------------------------------------------------

void GrackleChemistryData::set_int_(const std::string& field,
                                    const int& value) noexcept
{
  int* field_ptr = local_chemistry_data_access_int(ptr_.get(), field.c_str());
  if (!field_ptr){ field_err_("GrackleChemistryData::set_int_", field); }
  *field_ptr = value;
}

//----------------------------------------------------------------------

void GrackleChemistryData::set_dbl_(const std::string& field,
                                       const double& value) noexcept
{
  double* field_ptr = local_chemistry_data_access_double(ptr_.get(),
                                                         field.c_str());
  if (!field_ptr){ field_err_("GrackleChemistryData::set_dbl_", field); }
  *field_ptr = value;
}

//----------------------------------------------------------------------

void GrackleChemistryData::set_str_(const std::string& field,
                                    const std::string& value) noexcept
{
  const char** field_ptr = local_chemistry_data_access_string(ptr_.get(),
                                                              field.c_str());
  if (!field_ptr){ field_err_("GrackleChemistryData::set_str_", field); }

  // deallocate the existing value (if applicable)
  if ((*field_ptr) != nullptr){
    // we are explicitly checking for a matching pointer address
    auto is_match = [field_ptr](std::unique_ptr<char[]> elem)
                    { return elem.get() == (*field_ptr); };
    auto rslt = std::find_if(str_allocs_.begin(), str_allocs_.end(), is_match);

    // rslt will ONLY be equal to str_allocs_.end() if *field_ptr refers to a
    // string literal (that was set as a default value)
    if (rslt != str_allocs_.end()){ str_allocs_.erase(rslt); }
  }

  // allocate a new c-string and copy data from value into it
  const std::size_t length = value.size() + 1;
  std::unique_ptr<char[]> new_alloc(new char[length]);
  std::memcpy(new_alloc.get(), value.c_str(), length);

  // update the field of the chemistry_data struct
  (*field_ptr) = new_alloc.get();

  // finally, move the newly allocated c-string into str_allocs_
  str_allocs_.push_back(std::move(new_alloc));
}


GrackleChemistryData static GrackleChemistryData::from_parameters
  (Parameters& p,
   std::unordered_map<std::string, int> ext_int_params,
   std::unordered_map<std::string, double> ext_dbl_params,
   std::unordered_map<std::string, std::string> ext_str_params,
   std::unordered_set<std::string> unrelated_params) noexcept
{



}
