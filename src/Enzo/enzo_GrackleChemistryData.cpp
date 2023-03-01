// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_GrackleChemistryData.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com) 
/// @date     Feb 6 2023
/// @brief    [\ref Enzo] Implementation of the GrackleChemistryData class

#include "cello.hpp"

#include "enzo.hpp"

#include <cstring> // std::memcpy
#include <iterator> // std::input_iterator_tag
#include <limits> // std::numeric_limits
#include <type_traits> // std::is_same
#include <utility> // std::move
#include <functional> // std::function

//----------------------------------------------------------------------

#ifndef CONFIG_USE_GRACKLE
// provide a dummy definition of chemistry data
extern "C" { struct chemistry_data { int use_grackle; }; }
#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE

namespace{ // things within anonymous namespace aren't visible beyond this file

template<typename T>
std::string param_type_name(){
  if (std::is_same<T, int>::value) return "int";
  if (std::is_same<T, double>::value) return "double";
  if (std::is_same<T, std::string>::value) return "string";
  ERROR("param_type_name", "invalid template type");
}

//----------------------------------------------------------------------

/// returns the number of parameters of a given type. This is currently very
/// inefficient, but a function is being added to grackle to provide this exact
/// information in constant time
template<typename T>
int num_params(){
  std::function<const char*(std::size_t)> fn;
  if (std::is_same<T, int>::value) {
    fn = &param_name_int;
  } else if (std::is_same<T, double>::value) {
    fn = &param_name_double;
  } else if (std::is_same<T, std::string>::value) {
    fn = &param_name_string;
  } else {
    ERROR("num_params", "invalid template type");
  }

  std::size_t i = 0;
  while (fn(i) != nullptr) { i++; }

  ASSERT("num_params()", "unexpected result!", // this shouldn't happen!
         (std::size_t)(std::numeric_limits<int>::max()) >= i);

  return (int)i;
}

//----------------------------------------------------------------------

/// return ith parameter of the template type
template<typename T>
const char* param_name(int i){
  ASSERT("param_name", "arg is negative", i >= 0);
  if (std::is_same<T, int>::value)              {return param_name_int(i);}
  else if (std::is_same<T, double>::value)      {return param_name_double(i);}
  else if (std::is_same<T, std::string>::value) {return param_name_string(i);}
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
    std::string name_;
  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = std::string;
    using difference_type = long;
    using pointer = const std::string*;
    using reference = const std::string&;

    explicit iterator(long num) :
      num_(num), name_()
    {
      const char* tmp = param_name<T>(num);
      name_ = (tmp == nullptr) ? "" : tmp;
    }
    iterator& operator++() { *this = iterator(num_+1); return *this;}
    iterator operator++(int) { iterator out = *this; ++(*this); return out; }
    bool operator==(iterator other) const { return num_ == other.num_; }
    bool operator!=(iterator other) const { return !(other == *this); }
    reference operator*() const { return name_; } // should never be nullptr
  };

  iterator begin() {return iterator(0);}
  iterator end() { return iterator(num_params<T>()); }
};

//----------------------------------------------------------------------

/// if parameter is stored within chemistry_data, this returns the name of its
/// associated type. Otherwise, this returns an empty string.
///
/// This is slow and inefficient, but that's ok since it's intended for use
/// when formatting error messages
std::string type_name_of_param(const std::string& parameter){
  for (const std::string& name : GrackleParamNameRange<int>()){
    if (name == parameter){ return "int"; }
  }
  for (const std::string& name : GrackleParamNameRange<double>()){
    if (name == parameter){ return "double"; }
  }
  for (const std::string& name : GrackleParamNameRange<std::string>()){
    if (name == parameter){ return "string"; }
  }
  return "";
}

} // anonymous namespace

#endif /* CONFIG_USE_GRACKLE */

//----------------------------------------------------------------------

void GrackleChemistryData::param_err_(const std::string& func_name,
                                      const std::string& user_type,
                                      const std::string& field)
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::param_err_", "Grackle isn't linked");
#else

  std::string actual_type = type_name_of_param(field);

  if (actual_type != ""){
    ERROR3(func_name.c_str(),
           "the \"%s\" parameter holds a \"%s\" not a \"%s\"",
           field.c_str(), actual_type.c_str(), user_type.c_str());
  } else {
    ERROR1(func_name.c_str(), "No known parameter called \"%s\"",
           field.c_str());
  }

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

#ifndef CONFIG_USE_GRACKLE
GrackleChemistryData::GrackleChemistryData()
{ /* explicitly allow default construction */ }

#else

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
  // are both equivalent). For proof that this is the correct function to use,
  // see how the implementation of grackle's python extension
  //
  // there is an open pr that will let us replace this with a new function
  // called local_initialize_chemistry_parameters
  (*ptr_) = _set_default_chemistry_parameters();
}
#endif /* CONFIG_USE_GRACKLE */

//----------------------------------------------------------------------

GrackleChemistryData::~GrackleChemistryData() {}

//----------------------------------------------------------------------

GrackleChemistryData& GrackleChemistryData::operator=
(const GrackleChemistryData& other)
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::operator=", "Grackle isn't linked");
#else

  // first perform a direct assignment (this covers most fields correctly)
  (*ptr_) = *(other.ptr_);

  // Now, we manually copy over the strings. We explicitly use set_string to
  // ensure that lifetimes of strings get managed correctly
  for (const std::string& name : GrackleParamNameRange<std::string>()){
    this->set<std::string>(name, other.get<std::string>(name));
  }

  return *this;

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

namespace{ // things within anonymous namespace aren't visible beyond this file

std::string value_repr_(int val) { return std::to_string(val); }
std::string value_repr_(double val) { return std::to_string(val); }
std::string value_repr_(std::string val) { return "\"" + val + "\""; }

template<class T>
void pup_grackle_typed_params_(PUP::er &p, GrackleChemistryData& data)
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("pup_grackle_typed_params_", "Grackle isn't linked");
#else
  const bool up = p.isUnpacking();

  // first a sanity check:
  std::string param_type = param_type_name<T>();
  p | param_type;

  if (up) {
    std::string expected_type = param_type_name<T>();
    ASSERT2("pup_grackle_typed_params_",
            "This function expects to unpack values of type \"%s\", but the "
            "next type to be unpacked actually has type \"%s\"",
            expected_type.c_str(), param_type.c_str(),
            expected_type == param_type);
  }

  // now let's pup the number of parameters of the current type
  int expected_n_params = num_params<T>();
  int n_params = expected_n_params;
  p | n_params;

  // when expected_n_params < n_params: we will abort with a error when we try
  // to unpack an unknown parameter down below
  if (up && (expected_n_params > n_params)) {
    WARNING3("pup_grackle_typed_params_",
             "The function expects to unpack up to %d %s parameters, but only "
             "%d such parameters can be unpacked. Thus, the unpacked object "
             "may produce slightly different behavior from the original",
             expected_n_params, param_type.c_str(), n_params);
  }

  if (up){
    for (int i = 0; i < n_params; i++){
      std::string parameter;
      p | parameter;
      T value;
      p | value;

      if (!data.try_set<T>(parameter, value)){
        std::string err_prefix = ("Unpacked parameter \"" + parameter +
                                  "\" with an associated " + param_type +
                                  " value of " + value_repr_(value));
        std::string actual_type = type_name_of_param(parameter);
        std::string err_suffix = (actual_type == "") ?
          "isn't known" : ("should have a value of type " + actual_type);
        ERROR2("pup_grackle_typed_params_", "%s - this parameter %s",
               err_prefix.c_str(), err_suffix.c_str());
      }
    }
  } else {
    for (std::string parameter : GrackleParamNameRange<T>()){
      p | parameter;
      T value = data.get<T>(parameter);
      p | value;
    }
  }

#endif /* CONFIG_USE_GRACKLE */
}

} // anonymous namespace

//----------------------------------------------------------------------

void GrackleChemistryData::pup (PUP::er &p){
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::pup", "Grackle isn't linked");
#else
  // this function assumes that when this function is used for unpacking that
  // GrackleChemistryData's default constructor was used to initialize the
  // default values stored in this. This is EXTREMELY important when unpacking
  // data originally was serialized by a version of Enzo-E that was linked
  // against a different version of Grackle.
  
  // in the future, we may want to pup version information about Grackle (held
  // in the struct returned by get_grackle_version()). Then, during unpacking,
  // we could raise an:
  // - advisory warning informing the user the current version number doesn't
  //   exactly match the unpacked version information
  // - error if the current version number is smaller than the unpacked version
  //   number

  // Note: we could definitely speed this up when p.isMigration() evaluates to
  // true (if we know that the grackle library version will be identical on all
  // PEs, we could reduce message size by not pupping parameter names)
  pup_grackle_typed_params_<int>(p, *this);
  pup_grackle_typed_params_<double>(p, *this);
  pup_grackle_typed_params_<std::string>(p, *this);

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

std::pair<int,bool> GrackleChemistryData::try_get_int_
(const std::string& field) const noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::get_int_", "Grackle isn't linked");
#else

  int* field_ptr = local_chemistry_data_access_int
    (const_cast<chemistry_data*>(ptr_.get()), field.c_str());;
  if (!field_ptr) { return {0, false}; }
  return {*field_ptr, true};

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

std::pair<double,bool> GrackleChemistryData::try_get_dbl_
(const std::string& field) const noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::get_dbl_", "Grackle isn't linked");
#else

  double* field_ptr = local_chemistry_data_access_double
    (const_cast<chemistry_data*>(ptr_.get()), field.c_str());
  if (!field_ptr) { return {0.0, false}; }
  return {*field_ptr, true};

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

std::pair<std::string,bool> GrackleChemistryData::try_get_str_
(const std::string& field) const noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::get_str_", "Grackle isn't linked");
#else

  // we declare field_ptr as a double pointer to const characters to indicate
  // that the underlying values shouldn't be directly mutated through this
  // pointer (since the underlying data may be a string literal). For reasons
  // highlighted at:   https://c-faq.com/ansi/constmismatch.html
  // we must declare it as `const char * const *` instead of `const char **`
  const char * const * field_ptr = local_chemistry_data_access_string
    (const_cast<chemistry_data*>(ptr_.get()), field.c_str());
  if (!field_ptr) { return {"", false}; }
  return { ((*field_ptr) == nullptr) ? "" : *field_ptr,
           true };

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

bool GrackleChemistryData::try_set_int_(const std::string& field,
                                        const int& value) noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::try_set_int_", "Grackle isn't linked");
#else

  int* field_ptr = local_chemistry_data_access_int(ptr_.get(), field.c_str());
  if (!field_ptr) { return false; }
  *field_ptr = value;
  return true;

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

bool GrackleChemistryData::try_set_dbl_(const std::string& field,
                                        const double& value) noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::try_set_dbl_", "Grackle isn't linked");
#else

  double* field_ptr = local_chemistry_data_access_double(ptr_.get(),
                                                         field.c_str());
  if (!field_ptr) { return false; }
  *field_ptr = value;
  return true;

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

bool GrackleChemistryData::try_set_str_(const std::string& field,
                                        const std::string& value) noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::try_set_str_", "Grackle isn't linked");
#else

  // NOTE: we should NOT directly modify characters held by field_ptr
  char ** field_ptr = local_chemistry_data_access_string(ptr_.get(),
                                                         field.c_str());

  if (!field_ptr) { return false; }

  // deallocate the existing value (if applicable)
  if ((*field_ptr) != nullptr){

    // check whether *field_ptr matches the address of any pointers held within
    // str_allocs_, if so delete that pointer from str_allocs_
    for (std::size_t i = 0; i < str_allocs_.size(); i++){

      if (str_allocs_[i].get() == *field_ptr){
        // we will only enter this part of the loop if *field_ptr doesn't refer
        // to a string literal (that could have been set as a default value)
        str_allocs_.erase(str_allocs_.begin() + i);
        break;
      }
    }

  }

  // allocate a new c-string and copy data from value into it
  const std::size_t length = value.size() + 1;
  std::unique_ptr<char[]> new_alloc(new char[length]);
  std::memcpy(new_alloc.get(), value.c_str(), length);

  // update the field of the chemistry_data struct
  (*field_ptr) = new_alloc.get();

  // finally, move the newly allocated c-string into str_allocs_
  str_allocs_.push_back(std::move(new_alloc));
  return true;

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

GrackleChemistryData GrackleChemistryData::from_parameters
(Parameters& p, const str_vec_t& parameter_groups,
 const std::unordered_set<std::string>& forbid_leaf_names,
 const std::unordered_set<std::string>& ignore_leaf_names) noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleChemistryData::from_parameters", "Grackle isn't linked");
#else

  // save the groups for now...
  std::vector<std::string> group_v;
  for (int i = 0; i < p.group_count(); i++) { group_v.push_back(p.group(i)); }
  p.group_clear();

  std::string config_name_prefix = "";
  for (const std::string& part : parameter_groups){
    config_name_prefix.append(part + ":");
    p.group_push(part);
  }
  const str_vec_t leaf_parameter_names = p.leaf_parameter_names();
  p.group_clear();

  // restore the groups
  for (const std::string e : group_v) { p.group_push(e); }

  // ToDo: check that there aren't any groups with the shared prefix

  auto general_param_err =
    [config_name_prefix](const std::string& name, const std::string& user_type)
    {
      std::string grackle_param_type = type_name_of_param(name);
      const std::string full_conf_name = config_name_prefix + name;

      if (grackle_param_type == "") {
        ERROR2("GrackleChemistryData::from_parameters",
               "Invalid parameter specified in config file: %s. The linked "
               "Grackle library does NOT have a corresponding parameter: %s",
               full_conf_name.c_str(), name.c_str());
      } else {
        std::string tmp = (user_type == "") ? "incompatibly" : user_type;
        ERROR4("GrackleChemistryData::from_parameters",
               "The linked Grackle library expects the %s parameter to have "
               "a value of %s type. It can NOT be assigned the %s-typed value "
               "specified for the %s parameter in the config file",
               name.c_str(), grackle_param_type.c_str(),
               tmp.c_str(), full_conf_name.c_str());
      }
    };

  GrackleChemistryData out;

  for (const std::string& name : leaf_parameter_names){
    std::string full_conf_name = config_name_prefix + name;

    if (forbid_leaf_names.find(name) != forbid_leaf_names.end()) {
      ERROR1("GrackleChemistryData::from_parameters",
             "The config file is forbidden from having a parameter called %s",
             full_conf_name.c_str());
    } else if (ignore_leaf_names.find(name) != ignore_leaf_names.end()) {
      continue;
    }

    switch (p.type(full_conf_name)){
      case parameter_logical: {
        bool value = p.value_logical(full_conf_name);
        if (out.try_set<int>(name, value)) {break;}
        general_param_err(name, "logical");
      }
      case parameter_integer: {
        int value = p.value_integer(full_conf_name);
        if (out.try_set<int>(name, value)) {break;}
        general_param_err(name, "integer");
      }
      case parameter_float: {
        double value = p.value_float(full_conf_name);
        if (out.try_set<double>(name, value)) {break;}
        general_param_err(name, "float");
      }
      case parameter_string: {
        std::string value = p.value_string(full_conf_name);
        if (out.try_set<std::string>(name,value)) {break;}
        general_param_err(name, "string");
      }
      default:
        general_param_err(name, "");
    }
  }

  return out;

#endif /* CONFIG_USE_GRACKLE */
}

