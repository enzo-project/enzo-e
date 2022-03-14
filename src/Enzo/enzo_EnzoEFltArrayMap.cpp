// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEFltArrayMap.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon Oct 26 2020
/// @brief    [\ref Enzo] Implementation of the EnzoEFltArrayMap class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

namespace{ // define some local helper functions

  std::map<std::string, unsigned int> init_str_index_map_
    (const std::vector<std::string> &keys)
  {
    std::map<std::string, unsigned int> out;
    for (std::size_t i = 0; i < keys.size(); i++){
      const std::string& key = keys[i];
      ASSERT1("init_str_index_map_",
              "There can't be more than 1 key called %s",
              key.c_str(), out.find(key) == out.end());
      out[key] = (unsigned int)i;
    }
    return out;
  }

  void confirm_shared_shape_(const char* func_name,
			     const std::vector<std::string> &keys,
                             const std::vector<EFlt3DArray> &arrays){
    if (arrays.size() > 0){
      int ref_mz = arrays[0].shape(0);
      int ref_my = arrays[0].shape(1);
      int ref_mx = arrays[0].shape(2);
      for (std::size_t i = 1; i < arrays.size(); i++){
        int mz = arrays[i].shape(0);
        int my = arrays[i].shape(1);
        int mx = arrays[i].shape(2);

        ASSERT8(func_name,
                ("The shapes, (mz, my, mx), of the %s and %s arrays are "
                 "(%d, %d, %d) and (%d, %d, %d). The shapes must be the same"),
                keys[0].c_str(), keys[i].c_str(),
                ref_mz, ref_my, ref_mx, mz, my, mx,
                ((ref_mz == mz) && (ref_my == my) && (ref_mx == mx)));
      }
    }
  }

}

//----------------------------------------------------------------------

EnzoEFltArrayMap::EnzoEFltArrayMap(std::string name,
                                   const std::vector<std::string> &keys,
                                   const std::array<int,3>& shape)
  : name_(name)
{
  // copy the contents of keys to keys_
  keys_ = keys;

  // initialize str_index_map_
  str_index_map_ = init_str_index_map_(keys);

  // initialize arrays_
  arrays_.reserve(keys.size());
  for (const std::string& key : keys){
    arrays_.push_back(EFlt3DArray(shape[0], shape[1], shape[2]));
  }
}

//----------------------------------------------------------------------

EnzoEFltArrayMap::EnzoEFltArrayMap(std::string name,
                                   const std::vector<std::string> &keys,
                                   const std::vector<EFlt3DArray> &arrays)
  : name_(name)
{
  ASSERT2("EnzoEFltArrayMap::EnzoEFltArrayMap",
          "keys and arrays have lengths %zu and %zu. They should be the same",
          (std::size_t)keys.size(), (std::size_t)arrays.size(),
          keys.size() == arrays.size());

  // copy the contents of keys to keys_
  keys_ = keys;

  // initialize str_index_map_
  str_index_map_ = init_str_index_map_(keys);

  // validate that each entry in arrays_ has the same shape.
  confirm_shared_shape_("EnzoEFltArrayMap::EnzoEFltArrayMap", keys, arrays);

  // copy the contents of arrays to arrays_
  arrays_ = arrays;
}

//----------------------------------------------------------------------

bool EnzoEFltArrayMap::validate_key_order(const std::vector<std::string> &ref,
					  bool raise_err,
                                          bool allow_smaller_ref) const noexcept
{

  std::string name_string =
      (name_ == "") ? "" : (std::string(", \"") + name_ + std::string("\","));

  for (std::size_t i =0; i < std::min(ref.size(),keys_.size()); i++){
    bool equal = ref[i] == keys_[i];
    if (!equal && raise_err){
      print_summary();
      ERROR4("EnzoEFltArrayMap::validate_key_order",
	     ("keys of ArrayMap%s don't match expectations. At index %d, the "
	      "key is \"%s\" when it's expected to be \"%s\"\n"),
	     name_string.c_str(), (int)i, keys_[i].c_str(), ref[i].c_str());
    } else if (!equal){
      return false;
    }
  }

  if ((!allow_smaller_ref) && (ref.size() != keys_.size())){
    if (raise_err){
      print_summary();
      ERROR3("EnzoEFltArrayMap::validate_key_order",
	     "ArrayMap%s doesn't have the expected number of keys. It has %d "
	     "keys. It's expected to have %d",
             name_string.c_str(), (int)ref.size(), (int)keys_.size());
    } else {
      return false;
    }
  }
  return true;
}

//----------------------------------------------------------------------

CelloArray<enzo_float, 3> EnzoEFltArrayMap::at_(const std::string& key)
  const noexcept
{
  auto result = str_index_map_.find(key);
  if (result == str_index_map_.cend()){
    ERROR1("EnzoEFltArrayMap::at_", "map doesn't contain the key: \"%s\"",
           key.c_str());
  }
  return arrays_[result->second];
}

//----------------------------------------------------------------------

CelloArray<enzo_float, 3> EnzoEFltArrayMap::at_(const std::size_t index)
  const noexcept
{
  ASSERT("EnzoEFltArrayMap::at_",
         "index must be less than or equal to the length of the ArrayMap.",
         index < size());
  return arrays_[index];
}

//----------------------------------------------------------------------

EFlt3DArray exclude_stale_cells_(const EFlt3DArray &arr, int stale_depth)
{
  ASSERT("exclude_stale_cells_","each dim of arr must exceed 2*stale_depth.",
	 2*stale_depth < arr.shape(0) && 2*stale_depth < arr.shape(1) &&
	 2*stale_depth < arr.shape(2));

  return arr.subarray(CSlice(stale_depth, -stale_depth),
		      CSlice(stale_depth, -stale_depth),
		      CSlice(stale_depth, -stale_depth));
}


CelloArray<enzo_float, 3> EnzoEFltArrayMap::get_(const std::string& key,
                                                 int stale_depth) const noexcept
{
  ASSERT("EnzoEFltArrayMap::get_", "stale_depth must be >= 0",
         stale_depth >= 0);
  CelloArray<enzo_float, 3> arr = this->at_(key);
  return (stale_depth > 0) ? exclude_stale_cells_(arr,stale_depth) : arr;
}

//----------------------------------------------------------------------

void EnzoEFltArrayMap::print_summary() const noexcept
{
  std::size_t my_size = size();
  if (my_size == 0){
    if (name_ == ""){
      CkPrintf("Nameless Empty Array Map\n");
    } else {
      CkPrintf("\"%s\" Empty Array Map\n", name_.c_str());
    }
    return;
  }

  if (name_ == ""){
    CkPrintf("Nameless Array Map");
  } else {
    CkPrintf("\"%s\" Array Map", name_.c_str());
  }

  CkPrintf(": entry_shape = (%d, %d, %d)\n{",
           array_shape(0), array_shape(1), array_shape(2));

  for ( std::size_t i = 0; i < my_size; i++){
    if (i != 0){
      CkPrintf(",\n ");
    }

    EFlt3DArray array = arrays_[i];
    CkPrintf("\"%s\" : EFlt3DArray(%p)",
             keys_[i].c_str(), (void*)array.data());
  }
  CkPrintf("}\n");
  fflush(stdout);
}

//----------------------------------------------------------------------

int EnzoEFltArrayMap::array_shape(unsigned int dim) const noexcept{
  if (size() == 0){
    ERROR("EnzoEFltArrayMap::array_shape",
          "EnzoEFltArrayMap contains 0 arrays");
  }
  return arrays_[0].shape(dim);
}

//----------------------------------------------------------------------

// There might be some benefit to not directly constructing subarrays and lazily
// evaluating them instead in the output array
static inline std::vector<EFlt3DArray> make_subarrays_
(const std::vector<EFlt3DArray>& arrays,
 const CSlice &slc_z, const CSlice &slc_y, const CSlice &slc_x)
{
  std::vector<EFlt3DArray> out;
  out.reserve(arrays.size());
  for (const EFlt3DArray& arr : arrays){
    out.push_back(arr.subarray(slc_z, slc_y, slc_x));
  }
  return out;
}

EnzoEFltArrayMap EnzoEFltArrayMap::subarray_map(const CSlice &slc_z,
                                                const CSlice &slc_y,
                                                const CSlice &slc_x,
                                                const std::string& name)
{
  return EnzoEFltArrayMap(name, str_index_map_, keys_,
                          make_subarrays_(arrays_, slc_z, slc_y, slc_x));
}

const EnzoEFltArrayMap EnzoEFltArrayMap::subarray_map
(const CSlice &slc_z, const CSlice &slc_y, const CSlice &slc_x,
 const std::string& name) const
{
  return EnzoEFltArrayMap(name, str_index_map_, keys_,
                          make_subarrays_(arrays_, slc_z, slc_y, slc_x));
}

//----------------------------------------------------------------------

void EnzoEFltArrayMap::validate_invariants_() const noexcept
{
  ASSERT("EnzoEFltArrayMap::validate_invariants_",
         "str_index_map_, keys_, and arrays_ don't have the same length",
         keys_.size() == arrays_.size() &&
         keys_.size() == str_index_map_.size());

  for (std::size_t i = 0; i < keys_.size(); i++){
    auto result = str_index_map_.find(keys_[i]);
    ASSERT1("EnzoEFltArrayMap::validate_invariants_",
            "str_index_map_ is missing the key, \"%s\"", keys_[i].c_str(),
            result != str_index_map_.cend());

    ASSERT1("EnzoEFltArrayMap::validate_invariants_",
            ("str_index_map_ and keys_ disagree on the location of \"%s\" in "
             "the ordering"), keys_[i].c_str(),
            result->second == i);
  }
  
  // validate that each entry in arrays_ has the same shape.
  confirm_shared_shape_("EnzoEFltArrayMap::validate_invariants_",
                        keys_, arrays_);
}
