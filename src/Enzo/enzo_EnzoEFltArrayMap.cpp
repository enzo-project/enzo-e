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
  if (keys.size() > 0){
    int ref_mz = arrays[0].shape(0);
    int ref_my = arrays[0].shape(1);
    int ref_mx = arrays[0].shape(2);
    for (std::size_t i = 1; i < arrays.size(); i++){
      int mz = arrays[i].shape(0);
      int my = arrays[i].shape(1);
      int mx = arrays[i].shape(2);

      ASSERT8("EnzoEFltArrayMap::EnzoEFltArrayMap",
              ("The shapes, (mz, my, mx), of the %s and %s arrays are "
               "(%d, %d, %d) and (%d, %d, %d). The shapes must be the same"),
              keys[0].c_str(), keys[i].c_str(),
              ref_mz, ref_my, ref_mx, mz, my, mx,
              ((ref_mz == mz) && (ref_my == my) && (ref_mx == mx)));
    }
  }

  // copy the contents of arrays to arrays_
  arrays_ = arrays;
}

//----------------------------------------------------------------------

const EFlt3DArray& EnzoEFltArrayMap::at(const std::string& key) const noexcept
{
  auto result = str_index_map_.find(key);
  if (result == str_index_map_.cend()){
    ERROR1("EnzoEFltArrayMap::at", "map doesn't contain the key: \"%s\"",
           key.c_str());
  }
  return arrays_[result->second];
}

//----------------------------------------------------------------------

const EFlt3DArray& EnzoEFltArrayMap::at(const std::size_t index) const noexcept
{
  ASSERT("EnzoEFltArrayMap::at",
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


EFlt3DArray EnzoEFltArrayMap::get(const std::string& key,
                                  int stale_depth) const noexcept
{
  ASSERT("EnzoEFltArrayMap::get", "stale_depth must be >= 0",
         stale_depth >= 0);
  const EFlt3DArray& arr = this->at(key);
  if (stale_depth > 0){
    return exclude_stale_cells_(arr,stale_depth);
  } else {
    return arr;
  }
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
           arrays_[0].shape(0), arrays_[0].shape(1), arrays_[0].shape(2));

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
