// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEFltArrayMap.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon Oct 26 2020
/// @brief    [\ref Enzo] Implementation of the EnzoEFltArrayMap class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

const EFlt3DArray& EnzoEFltArrayMap::at(const std::string& key) const noexcept
{
  auto result = map_.find(key);
  if (result == map_.cend()){
    ERROR1("EnzoEFltArrayMap::at", "map doesn't contain the key: \"%s\"",
           key.c_str());
  }
  return result->second;
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
    CkPrintf("Nameless Array Map\n");
  } else {
    CkPrintf("\"%s\" Array Map\n", name_.c_str());
  }

  int i = 0;
  CkPrintf("{");

  for ( const auto &pair : map_ ) {
    if (i != 0){
      CkPrintf(",\n ");
    }
    CkPrintf("\"%s\" : EFlt3DArray(%p, %d, %d, %d), owners: %ld",
             pair.first.c_str(),
             (void*)pair.second.shared_data_.get(),
             (int)pair.second.shape(0),
             (int)pair.second.shape(1),
             (int)pair.second.shape(2),
             pair.second.shared_data_.use_count());
    i++;
  }
  CkPrintf("}\n");
  fflush(stdout);
}
