// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEFltArrayMap.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon Oct 26 2020
/// @brief    [\ref Enzo] Implementation of the EnzoEFltArrayMap class

#include "Enzo/utils/utils.hpp"
#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"

//----------------------------------------------------------------------

EnzoEFltArrayMap::EnzoEFltArrayMap(std::string name,
                                   const std::vector<std::string> &keys,
                                   const std::array<int,3>& shape)
  : name_(name),
    str_index_map_(keys),
    arrays_(keys.size(), shape)
{ }

//----------------------------------------------------------------------

EnzoEFltArrayMap::EnzoEFltArrayMap(std::string name,
                                   const std::vector<std::string> &keys,
                                   const std::vector<EFlt3DArray> &arrays)
  : name_(name),
    str_index_map_(keys),
    arrays_(arrays)
{
  ASSERT2("EnzoEFltArrayMap::EnzoEFltArrayMap",
          "keys and arrays have lengths %zu and %zu. They should be the same",
          (std::size_t)keys.size(), (std::size_t)arrays.size(),
          keys.size() == arrays.size());
}

//----------------------------------------------------------------------

bool EnzoEFltArrayMap::validate_key_order(const std::vector<std::string> &ref,
					  bool raise_err,
                                          bool allow_smaller_ref) const noexcept
{

  std::string name_string =
      (name_ == "") ? "" : (std::string(", \"") + name_ + std::string("\","));

  for (std::size_t i =0; i < std::min(ref.size(),str_index_map_.size()); i++){
    const std::string k = str_index_map_.key(i);
    bool equal = ref[i] == k;
    if (!equal && raise_err){
      print_summary();
      ERROR4("EnzoEFltArrayMap::validate_key_order",
	     ("keys of ArrayMap%s don't match expectations. At index %d, the "
	      "key is \"%s\" when it's expected to be \"%s\"\n"),
	     name_string.c_str(), (int)i, k.c_str(), ref[i].c_str());
    } else if (!equal){
      return false;
    }
  }

  if ((!allow_smaller_ref) && (ref.size() != str_index_map_.size())){
    if (raise_err){
      print_summary();
      ERROR3("EnzoEFltArrayMap::validate_key_order",
	     "ArrayMap%s doesn't have the expected number of keys. It has %d "
	     "keys. It's expected to have %d",
             name_string.c_str(), (int)ref.size(), (int)str_index_map_.size());
    } else {
      return false;
    }
  }
  return true;
}

//----------------------------------------------------------------------

CelloView<enzo_float, 3> EnzoEFltArrayMap::at_(const std::string& key)
  const noexcept
{ return arrays_[str_index_map_.at(key)]; }

//----------------------------------------------------------------------

CelloView<enzo_float, 3> EnzoEFltArrayMap::at_(const std::size_t index)
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


CelloView<enzo_float, 3> EnzoEFltArrayMap::get_(const std::string& key,
                                                 int stale_depth) const noexcept
{
  ASSERT("EnzoEFltArrayMap::get_", "stale_depth must be >= 0",
         stale_depth >= 0);
  CelloView<enzo_float, 3> arr = this->at_(key);
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
    std::string key = str_index_map_.key(i);
    CkPrintf("\"%s\" : EFlt3DArray(%p)", key.c_str(), (void*)array.data());
  }
  CkPrintf("}\n");
  fflush(stdout);
}

//----------------------------------------------------------------------

EnzoEFltArrayMap EnzoEFltArrayMap::subarray_map(const CSlice &slc_z,
                                                const CSlice &slc_y,
                                                const CSlice &slc_x,
                                                const std::string& name)
{
  return EnzoEFltArrayMap(name, str_index_map_,
                          arrays_.subarray_collec(slc_z, slc_y, slc_x));
}

const EnzoEFltArrayMap EnzoEFltArrayMap::subarray_map
(const CSlice &slc_z, const CSlice &slc_y, const CSlice &slc_x,
 const std::string& name) const
{
  return EnzoEFltArrayMap(name, str_index_map_,
                          arrays_.subarray_collec(slc_z, slc_y, slc_x));
}

//----------------------------------------------------------------------

void EnzoEFltArrayMap::validate_invariants_() const noexcept
{
  // several invariants are alread enforced:
  // - StringIndRdOnlyMap implicitly enforces that there aren't any duplicate
  //   keys, and that a unique key is associated with each integer index from 0
  //   through (str_index_map_.size() - 1)
  // - ViewCollec enforces that all arrays have the same shape
  ASSERT("EnzoEFltArrayMap::validate_invariants_",
         "str_index_map_ and arrays_ don't have the same length",
         arrays_.size() == str_index_map_.size());
}
