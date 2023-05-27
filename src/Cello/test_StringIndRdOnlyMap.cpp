// See LICENSE_CELLO file for license and copyright information

/// @file     test_StringIndRdOnlyMap.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-04-01
/// @brief    Test program for StringIndRdOnlyMap

#include "main.hpp"
#include "test.hpp"
#include "view.hpp"

//----------------------------------------------------------------------

bool consistent_key_order_(const std::vector<std::string>& ref_order,
                           const StringIndRdOnlyMap& string_ind_map)
{
  if (ref_order.size() != string_ind_map.size()){
    return false;
  }

  for (std::size_t i = 0; i < string_ind_map.size(); i++){
    const std::string& key = ref_order[i];
    
    // first, check that the key is even contained by string_ind_map
    if (!string_ind_map.contains(key)) { return false; }
    
    // next, check the value associated with key using the operator[] method
    if (string_ind_map[key] != i) { return false; }

    // then, check the value associated with key using the at method
    if (string_ind_map.at(key) != i) { return false; }

    // finally, check result with key method
    if (key != string_ind_map.key(i)) { return false; }
  }
  return true;
}

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("StringIndRdOnlyMap");

  std::vector<std::string> key_order_1 = {"density", "velocity_x", "velocity_y",
                                          "velocity_z", "total_energy"};
  std::vector<std::string> key_order_2 = {"total_energy", "density",
                                          "velocity_x", "velocity_y",
                                          "velocity_z"};
  std::vector<std::string> key_order_3 = {"bfield_x", "internal_energy"};

  // it would probably be better to isolate which methods are being tested at
  // a given time. This first test checks a bunch of things.
  StringIndRdOnlyMap str_ind_map_primary(key_order_1);

  unit_assert (consistent_key_order_(key_order_1, str_ind_map_primary));
  // this is effectively a sanity check
  unit_assert (!consistent_key_order_(key_order_2, str_ind_map_primary));
  unit_assert (!consistent_key_order_(key_order_3, str_ind_map_primary));


  // now check copy constructor:
  {
    StringIndRdOnlyMap str_ind_map_secondary = str_ind_map_primary;

    unit_assert (consistent_key_order_(key_order_1, str_ind_map_secondary));
    // check that the original is unaffected by the copy:
    unit_assert (consistent_key_order_(key_order_1, str_ind_map_primary));
  }

  // now check copy assignment:
  {
    StringIndRdOnlyMap str_ind_map_other(key_order_3);
    unit_assert (consistent_key_order_(key_order_3, str_ind_map_other));

    str_ind_map_other = str_ind_map_primary;
    unit_assert (!consistent_key_order_(key_order_3, str_ind_map_other));
    unit_assert (consistent_key_order_(key_order_1, str_ind_map_other));

    // check that the original is unaffected by the copy:
    unit_assert (consistent_key_order_(key_order_1, str_ind_map_primary));
  }

  //----------------------------------------------------------------------
  unit_finalize();
  //----------------------------------------------------------------------

  exit_();
}
PARALLEL_MAIN_END
