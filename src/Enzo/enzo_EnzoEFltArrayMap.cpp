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

EnzoEFltArrayMap EnzoEFltArrayMap::from_grouping
(Block * block, Grouping& grouping, const std::vector<std::string>& groups,
 int dim) noexcept
{
  EnzoEFltArrayMap out;

  EnzoFieldArrayFactory array_factory(block,0);

  for (const std::string& group_name : groups){
    int num_fields = grouping.size(group_name);
    for (int field_ind=0; field_ind<num_fields; field_ind++){
      std::string field_name = grouping.item(group_name,field_ind);

      if (out.contains(field_name)){
        ERROR1("EnzoEFltArrayMap::from_grouping",
               "EnzoEFltArrayMap can't hold more than one field called \"%s\"",
               field_name.c_str());
      }

      if (dim == -1){
        out[field_name] = array_factory.from_name(field_name);
      } else {
        out[field_name] = array_factory.assigned_center_from_name(field_name,
                                                                  dim);
      }
    }
  }
  return out;
}

//----------------------------------------------------------------------

void EnzoEFltArrayMap::print_summary() const noexcept
{
  std::size_t my_size = size();
  if (my_size == 0){
    CkPrintf("Empty array");
    return;
  }

  int i = 0;
  CkPrintf("{");
  
  for ( const auto &pair : map_ ) {
    if (i != 0){
      CkPrintf(",\n ");
    }
    CkPrintf("\"%s\" : EFlt3DArray(%p, %d, %d, %d)", pair.first.c_str(),
             (void*)pair.second.shared_data_.get(),
             (int)pair.second.shape(0),
             (int)pair.second.shape(1),
             (int)pair.second.shape(2));
    i++;
  }
  CkPrintf("}\n");
  fflush(stdout);
}
