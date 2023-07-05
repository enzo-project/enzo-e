// See LICENSE_CELLO file for license and copyright information

/// @file     test_ViewMap.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-07-03
/// @brief    Test program for the ViewMap

#include "main.hpp"
#include "test.hpp"
#include "view.hpp"

// instruct the ViewTestRoutines to use unit_assert
#define VIEW_TEST_ROUTINES_EMPLOY_UNIT_ASSERT 1
#include "test_ViewTestRoutines.hpp"

// VIEW_TEST_ASSERT was defined in test_ViewTestRoutines.hpp
// - we may want to move this definition to test_Unit.hpp (and change the name)
// - it prints out a detailed error message upon failure (just like ASSERT),
//   but also hooks into unit-testing

#include <array>
#include <utility> // for std::pair
#include <vector>

//----------------------------------------------------------------------

class SimpleInitializationTests{
  // This class holds tests used to check the simple initialization of
  // ViewMap. Unfortunately, it is not really possible to completely
  // disentangle the following set of features and test them completely
  // independently of each other:
  // 1. direct initialization of an instance that wraps existing CelloView
  // 2. direct initialization of an instance which manages its own data (and is
  //    contiguous)
  // 3. retrieval operations (we could pick any retrieval operation - in this
  //    context, we pick the at method)
  // 4. the name method

public:

  void default_constructor(){
    unit_func_quiet("ViewMap", "ViewMap()");
    ViewMap<int> nameless_view_map;
    VIEW_TEST_ASSERT("SimpleInitializationTests::default_constructor",
                     "dflt constructor should produce viewmap with empty name",
                     nameless_view_map.name() == "");
    VIEW_TEST_ASSERT("SimpleInitializationTests::default_constructor",
                     "dflt constructor should produce viewmap with 0 elements",
                     nameless_view_map.size() == 0);
  }

  void default_named_constructor(){
    unit_func_quiet("ViewMap", "ViewMap(std::string)");
    ViewMap<int> view_map("my_ViewMap");
    VIEW_TEST_ASSERT("SimpleInitializationTests::default_named_constructor",
                     "default named constructor set an invalid name",
                     view_map.name() == "my_ViewMap");
    VIEW_TEST_ASSERT("SimpleInitializationTests::default_named_constructor",
                     "default named constructor should produce viewmap with 0 "
                     "elements",
                     view_map.size() == 0);
  }

  template<typename T>
  void require_all_zero_(ViewMap<T>& view_map,
                         const std::string& expected_name,
                         const std::vector<std::string>& keys,
                         const std::array<int,3>& shape)
  {
    // the purpose of this method is to do comparisons while using as few
    // methods of ViewMap as possible since those other methods haven't been
    // validated separately
    //
    // (we are free to use any methods of CelloView)
    for (const std::string& key: keys){
      CelloView<int,3> view = view_map.at(key);
      assert_shape3D(view, shape);

      bool all_equal = true;
      for (int iz = 0; iz < shape[0]; iz++) {
        for (int iy = 0; iy < shape[1]; iy++) {
          for (int ix = 0; ix < shape[2]; ix++) {
            all_equal &= ( view(iz,iy,ix) == 0);
          }
        }
      }

      std::string msg = "The elements of key " + key + "are not all 0";

      VIEW_TEST_ASSERT("SimpleInitializationTests::_require_all_zero",
                       msg.c_str(), all_equal);
    }
  }

  void test_init_simple(bool contiguous){

    // all values in a CelloView instance are assumed to be zero-initialized

    std::vector<std::string> names = {"", "my_ViewMap"};
    std::vector<std::string> keys = {"density", "velocity_x", "total_energy"};
    std::array<int,3> shape = {4, 6, 8};

    for (auto name : names) {
      ViewMap<int> *view_map;
      std::vector<CelloView<int, 3>> view_vec; // only used in 1 case
      if (contiguous) {
        if (name == "") {
          unit_func_quiet("ViewMap", "ViewMap(std::vector<std::string>, std::array<int,3>)");

          view_map = new ViewMap<int>(keys, shape);
          this->require_all_zero_(*view_map, "", keys, shape);
        } else {
          unit_func_quiet("ViewMap", "ViewMap(std::string, std::vector<std::string>, std::array<int,3>)");

          view_map = new ViewMap<int>(name, keys, shape);
          this->require_all_zero_(*view_map, name, keys, shape);
        }

      } else {
        view_vec.push_back(CelloView<int, 3>(shape[0], shape[1], shape[2]));
        view_vec.push_back(CelloView<int, 3>(shape[0], shape[1], shape[2]));
        view_vec.push_back(CelloView<int, 3>(shape[0], shape[1], shape[2]));

        if (name == "") {
          unit_func_quiet("ViewMap", "ViewMap(std::vector<std::string>, std::vector<CelloView<T,3>>)");
          
          view_map = new ViewMap<int>(keys, view_vec);
          this->require_all_zero_(*view_map, "", keys, shape);
        } else {
          unit_func_quiet("ViewMap", "ViewMap(std::string, std::vector<std::string>, std::vector<CelloView<T,3>>)");
          
          view_map = new ViewMap<int>(name, keys, view_vec);
          this->require_all_zero_(*view_map, name, keys, shape);
        }
      }

      unit_assert(view_map->contiguous_arrays() == contiguous);
      delete view_map;
    }
  }

  void test_init_nonzero_wrapped(){
    unit_func_quiet("ViewMap", "ViewMap(std::vector<std::string>, std::vector<CelloView<T,3>>)");
    std::vector<std::string> keys = {"density", "velocity_x", "total_energy"};
    std::array<int,3> shape = {4, 6, 8};

    std::vector<CelloView<int, 3>> view_vec =
      { range_3Darray(shape[0], shape[1], shape[2], ValRange<int>(0, 1)),
        range_3Darray(shape[0], shape[1], shape[2], ValRange<int>(10000, 1)),
        range_3Darray(shape[0], shape[1], shape[2], ValRange<int>(20000, 1)) };

    ViewMap<int> view_map = ViewMap<int>(keys, view_vec);

    assert_shape3D(view_map.at("density"), shape);
    assert_allequal3D(view_map.at("density"), ValRange<int>(0, 1));

    assert_shape3D(view_map.at("velocity_x"), shape);
    assert_allequal3D(view_map.at("velocity_x"), ValRange<int>(10000, 1));

    assert_shape3D(view_map.at("total_energy"), shape);
    assert_allequal3D(view_map.at("total_energy"), ValRange<int>(20000, 1));

    // now let's make sure we reflect changes in the underlying array
    // (without modifying the array returned by the at method - we test that
    // separately)
    assign_range_3Darray(view_vec[1], ValRange<int>(17, -1));
    assert_allequal3D(view_map.at("velocity_x"), ValRange<int>(17, -1));
  }

  void run_tests(){
    default_constructor();
    default_named_constructor();
    test_init_simple(true);
    test_init_simple(false);
    // it is not possible to initialize non-zero managed memory
    test_init_nonzero_wrapped();
  }

};

//----------------------------------------------------------------------

template<typename T, bool IsContiguous>
class ManagedViewMap{
  // This template class is used to build and a manage the lifetime of a
  // ViewMap instance
public:

  ManagedViewMap(const std::vector<std::string>& keys,
                 const std::array<int,3>& shape)
  {
    if (IsContiguous) {
      view_map_ptr_ = new ViewMap<T>(keys, shape);
    } else {
      for (std::size_t i = 0; i < keys.size(); i++) {
        view_vec_.push_back(CelloView<T,3>(shape[0], shape[1], shape[2]));
      }
      view_map_ptr_ = new ViewMap<T>(keys, view_vec_);
    }
  }

  ~ManagedViewMap() { delete view_map_ptr_; }

  ManagedViewMap() = delete;
  ManagedViewMap(const ManagedViewMap&) = delete;
  ManagedViewMap& operator=(const ManagedViewMap&) = delete;

  // strictly speaking, move constructor/assignment doesn't NEED to be deleted
  ManagedViewMap(ManagedViewMap&&) = delete;
  ManagedViewMap& operator=(ManagedViewMap&&) = delete;

  ViewMap<T>& get() const {return *view_map_ptr_;}
  ViewMap<T>* get_ptr() const {return view_map_ptr_; }
  bool is_contiguous() const {return IsContiguous; }

  const std::vector<CelloView<T,3>>& view_vec() const {return view_vec_;}

private:

  // the following is always not a nullptr
  ViewMap<T>* view_map_ptr_;

  // the following is only conditionally used
  std::vector<CelloView<T,3>> view_vec_;
};

//----------------------------------------------------------------------

class ViewAccessTests{
  // testing functionality associated with accessing contained Views
public:

  template<bool Contiguous>
  void test_view_accessor(bool bracket_index){

    std::string func_name =
      (bracket_index) ? "operator[](int)" : "operator[](std::string)";
    Unit::instance()->set_func("ViewMap", func_name.c_str());

    const std::vector<std::string> keys =
      {"density", "velocity_x", "total_energy"};
    std::array<int,3> shape = {4, 6, 8};

    ManagedViewMap<int, Contiguous> manager(keys, shape);
    ViewMap<int>* view_map = manager.get_ptr();

    auto accessor = [=](std::string key) -> CelloView<int,3>
      {
        if (bracket_index){
          int index;
          if (key == "density") { index = 0; }
          else if (key == "velocity_x") { index = 1; }
          else if (key == "total_energy") { index = 2; }
          else {
            ERROR("ViewAccessTests::test_view_accessor", "invalid key");
          }
          return view_map->operator[](index);
        } else {
          return view_map->operator[](key);
        }
      };

    // 1st check that updates applied to CelloView instances retrieved with
    // ViewMap::at are seen when the CelloView is accessed by target accessor
    {
      CelloView<int, 3> dens_view = accessor("density");
      CelloView<int, 3> vx_view = accessor("velocity_x");
      CelloView<int, 3> etot_view = accessor("total_energy");
    
      assign_range_3Darray(view_map->at("density"), ValRange<int>(0, 1));
      assign_range_3Darray(view_map->at("velocity_x"),
                           ValRange<int>(10000, 1));
      assign_range_3Darray(view_map->at("total_energy"),
                           ValRange<int>(20000, 1));

      assert_allequal3D(dens_view, ValRange<int>(0, 1));
      assert_allequal3D(vx_view, ValRange<int>(10000, 1));
      assert_allequal3D(etot_view, ValRange<int>(20000, 1));
    }

    // next check that updates applied to CelloView instances retrieved by
    // target accessor are seen when the CelloView is accessed with ViewMap::at
    {
      CelloView<int, 3> dens_view = accessor("density");
      CelloView<int, 3> vx_view = accessor("velocity_x");
      CelloView<int, 3> etot_view = accessor("total_energy");
    
      assign_range_3Darray(dens_view, ValRange<int>(-100, -1));
      assign_range_3Darray(vx_view, ValRange<int>(15, 2));
      assign_range_3Darray(etot_view, ValRange<int>(-34, 7));

      assert_allequal3D(view_map->at("density"), ValRange<int>(-100, -1));
      assert_allequal3D(view_map->at("velocity_x"), ValRange<int>(15, 2));
      assert_allequal3D(view_map->at("total_energy"), ValRange<int>(-34, 7));
    }
  }

  void run_tests(){
    test_view_accessor<true>(true);
    test_view_accessor<true>(false);
    test_view_accessor<false>(true);
    test_view_accessor<false>(false);
  }
};

//----------------------------------------------------------------------

template <typename T>
bool consistent_key_order_(const std::vector<std::string>& ref_order,
                           const ViewMap<T>& view_map)
{
  if (ref_order.size() != view_map.size()){
    return false;
  }

  for (std::size_t i = 0; i < view_map.size(); i++){
    const std::string& key = ref_order[i];
    
    // first, check that the key is even contained by view_map
    if (!view_map.contains(key)) { return false; }
    
    // next, check the value associated with key using the operator[] method
    if (!view_map[key].is_alias(view_map[i])) { return false; }

    // then, check the value associated with key using the at method
    if (!view_map.at(key).is_alias(view_map[i])) { return false; }

  }
  return true;
}

//----------------------------------------------------------------------

class MiscOperations{
  // these are tests of misc methods.

public:

  void test_backing_array() {
    Unit::instance()->set_func("ViewMap", "get_backing_array");

    // we can't test this unless we are using a contiguous ViewMap
    
    const std::vector<std::string> keys = {"A", "B", "C"};
    const std::array<int,3> shape = {4, 6, 8};

    ViewMap<int> view_map(keys, shape);

    view_map.at("B")(2,3,4) = 3;

    // check that initial values of backing_array are what we expect
    CelloView<int, 4> backing_array = view_map.get_backing_array();
    {
      bool all_equal = true;
      for (int ikey = 0; ikey < 3; ikey++) {
        for (int iz = 0; iz < shape[0]; iz++) {
          for (int iy = 0; iy < shape[1]; iy++) {
            for (int ix = 0; ix < shape[2]; ix++) {
              int expected_val = 0;
              if ((ikey == 1) & (iz == 2) & (iy == 3) & (ix == 4)){
                expected_val = 3;
              }
              all_equal &= (backing_array(ikey, iz, iy, ix) == expected_val);
            }
          }
        }
      }

      unit_assert(all_equal);
    }


    // check that mutations to the backing array are reflected in view_map
    backing_array(0,1,2,3) = -52;
    backing_array(2,0,0,0) = 1;
    {
      bool all_equal = true;
      for (const std::string& key: keys) {
        CelloView<int,3> cur_view = view_map.at(key);

        for (int iz = 0; iz < shape[0]; iz++) {
          for (int iy = 0; iy < shape[1]; iy++) {
            for (int ix = 0; ix < shape[2]; ix++) {
              int expected_val = 0;
              if ((key == "A") & (iz == 1) & (iy == 2) & (ix == 3)) {
                expected_val = -52;
              } else if ((key == "B") & (iz == 2) & (iy == 3) & (ix == 4)){
                expected_val = 3;
              } else if ((key == "C") & (iz == 0) & (iy == 0) & (ix == 0)){
                expected_val = 1;
              }
              all_equal &= (cur_view(iz, iy, ix) == expected_val);
            }
          }
        }
      }

      unit_assert(all_equal);
      //ASSERT("MiscOperations::test_backing_array",
      //       "view_map has unexpected value after mutating backing_array",
      //       all_equal);
    }
  }

  template<bool isContiguous>
  void test_array_shape(){
    Unit::instance()->set_func("ViewMap", "array_shape");

    const std::vector<std::string> keys = {"A", "B", "C"};
    const std::array<int,3> shape = {4, 6, 8};

    ManagedViewMap<int, isContiguous> manager(keys, shape);
    ViewMap<int>& view_map = manager.get();

    unit_assert(view_map.array_shape(0) == 4);
    unit_assert(view_map.array_shape(1) == 6);
    unit_assert(view_map.array_shape(2) == 8);
  }

  void test_default_constructed_size() {
    Unit::instance()->set_func("ViewMap", "size");

    ViewMap<int> view_map;
    unit_assert(view_map.size() == 0);
  }

  template<bool isContiguous>
  void test_size() {
    Unit::instance()->set_func("ViewMap", "size");

    const std::array<int,3> shape = {4, 6, 8};

    {
      ManagedViewMap<int, isContiguous> manager({"my_key"}, shape);
      unit_assert(manager.get().size() == 1);
    }

    {
      ManagedViewMap<int, isContiguous> manager({"A", "B", "C"}, shape);
      unit_assert(manager.get().size() == 3);
    }
  }

  template<bool isContiguous>
  void test_misc_map_methods() {
    Unit::instance()->set_class("ViewMap");

    const std::vector<std::string> keys = {"A", "B", "C"};
    const std::array<int,3> shape = {4, 6, 8};

    ManagedViewMap<int, isContiguous> manager(keys, shape);
    unit_assert(consistent_key_order_(keys, manager.get()));
  }

  void run_tests(){
    test_backing_array();
    test_array_shape<true>();
    test_array_shape<false>();

    test_default_constructed_size();
    test_size<true>();
    test_size<false>();
    
    test_misc_map_methods<true>();
    test_misc_map_methods<false>();
  }
};

//----------------------------------------------------------------------


class CopyTests{
  // these are tests that check that check the assignments of ViewMaps to
  // variables.

private:

  template<typename ManagedViewMapType>
  void common_checks_(ViewMap<int> &dest_view_map,
                      ManagedViewMapType& src_manager,
                      const std::vector<std::string>& src_keys) {

    ViewMap<int>& src_view_map = src_manager.get();

    // perform a check related to map-methods
    unit_assert(consistent_key_order_(src_keys, dest_view_map));

    // now let's perform some checks to ensure that the shallow copy operations
    // work as expected

    assign_range_3Darray(src_view_map.at("B"), ValRange<int>(17, -1));
    assert_allequal3D(dest_view_map.at("B"), ValRange<int>(17, -1));

    assign_range_3Darray(dest_view_map.at("C"), ValRange<int>(0, 1));
    assert_allequal3D(src_view_map.at("C"), ValRange<int>(0, 1));

    if (src_manager.is_contiguous()) {
      // we could add some tests related to the backing array
    } else {
      assert_allequal3D(src_manager.view_vec()[0], 0);
      assert_allequal3D(src_manager.view_vec()[1], ValRange<int>(17, -1));
      assert_allequal3D(src_manager.view_vec()[2], ValRange<int>(0, 1));
    }
  }

public:

  template<bool Contiguous>
  void test_constructor() {
    Unit::instance()->set_func("ViewMap", "ViewMap(const ViewMap&)");

    // setup the src (that will be assigned)
    const std::vector<std::string> src_keys = {"A", "B", "C"};
    const std::array<int,3> src_shape = {4, 6, 8};
    ManagedViewMap<int, Contiguous> src_manager(src_keys, src_shape);

    // execute the copy Constructor
    ViewMap<int> dest_view_map(src_manager.get());

    // perform the checks:
    common_checks_(dest_view_map, src_manager, src_keys);
  }

  template<bool srcContiguous, bool destContiguous>
  void test_assignment() {
    Unit::instance()->set_func("ViewMap", "operator=(const ViewMap&)");

    // setup the src (that will be assigned)
    const std::vector<std::string> src_keys = {"A", "B", "C"};
    const std::array<int,3> src_shape = {4, 6, 8};
    ManagedViewMap<int, srcContiguous> src_manager(src_keys, src_shape);

    // setup the dest (that will be overwritten by the assignment)
    const std::vector<std::string> dest_init_keys = {"C", "D", "E"};
    const std::array<int,3> dest_init_shape = {3, 2, 1};
    ManagedViewMap<int, destContiguous> dest_manager(dest_init_keys,
                                                    dest_init_shape);
    ViewMap<int>& dest_view_map = dest_manager.get();

    // perform the actual assignment
    dest_view_map = src_manager.get();

    // perform the checks:
    common_checks_(dest_view_map, src_manager, src_keys);
  }

  void run_tests()
  {
    test_constructor<true>();
    test_constructor<false>();

    test_assignment<true, true>();
    test_assignment<true, false>();
    test_assignment<false, true>();
    test_assignment<false, false>();
  }
};

//----------------------------------------------------------------------

class MoveTests{
  // tests that check ViewMap's move constructor and move assignment

private:

  // the following gets invoked after the move constructor/assignment.
  //
  // If view_vec is empty, then we assume that the source view_map allocated
  // its own memory
  void common_checks_(ViewMap<int> &view_map,
                      const std::vector<std::string>& src_keys,
                      const std::vector<ValRange<int>>& expected_vranges,
                      const std::vector<CelloView<int,3>>& view_vec) {

    // perform a check related to map-methods
    unit_assert(consistent_key_order_(src_keys, view_map));

    for (std::size_t i = 0; i < expected_vranges.size(); i++){
      assert_allequal3D(view_map[i], expected_vranges[i]);

      if (view_vec.size() != 0) {
        assert_allequal3D(view_vec[i], expected_vranges[i]);
      } else {
        // it may make sense to test the backing array
      }
    }
  }

public:

  template<bool Contiguous>
  void test_constructor() {
    Unit::instance()->set_func("ViewMap", "ViewMap(ViewMap&&)");

    // setup the src (that will be assigned)
    const std::vector<std::string> src_keys = {"A", "B", "C"};
    const std::array<int,3> src_shape = {4, 6, 8};
    ManagedViewMap<int, Contiguous> src_manager(src_keys, src_shape);

    assign_range_3Darray(src_manager.get().at("B"), ValRange<int>(17, -1));

    // execute the move Constructor
    ViewMap<int> dest_view_map(std::move(src_manager.get()));

    assign_range_3Darray(dest_view_map.at("C"), ValRange<int>(0, 1));

    // perform the checks:
    common_checks_
      (dest_view_map, src_keys,
       {ValRange<int>(0,0), ValRange<int>(17,-1), ValRange<int>(0,1)},
       src_manager.view_vec());
  }


  template<bool srcContiguous, bool destContiguous>
  void test_assignment() {
    Unit::instance()->set_func("ViewMap", "operator=(ViewMap&&)");

    // setup the src (that will be assigned)
    const std::vector<std::string> src_keys = {"A", "B", "C"};
    const std::array<int,3> src_shape = {4, 6, 8};
    ManagedViewMap<int, srcContiguous> src_manager(src_keys, src_shape);

    // setup the dest (that will be overwritten by the assignment)
    const std::vector<std::string> dest_init_keys = {"C", "D", "E"};
    const std::array<int,3> dest_init_shape = {3, 2, 1};
    ManagedViewMap<int, destContiguous> dest_manager(dest_init_keys,
                                                    dest_init_shape);
    ViewMap<int>& dest_view_map = dest_manager.get();

    assign_range_3Darray(src_manager.get().at("B"), ValRange<int>(17, -1));

    // perform the actual assignment
    dest_view_map = std::move(src_manager.get());

    assign_range_3Darray(dest_view_map.at("C"), ValRange<int>(0, 1));

    // perform the checks:
    common_checks_
      (dest_view_map, src_keys,
       {ValRange<int>(0,0), ValRange<int>(17,-1), ValRange<int>(0,1)},
       src_manager.view_vec());
  }


  void run_tests()
  {
    test_constructor<true>();
    test_constructor<false>();

    test_assignment<true, true>();
    test_assignment<true, false>();
    test_assignment<false, true>();
    test_assignment<false, false>();
  }
};

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{
  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ViewMap");

  SimpleInitializationTests init_tests;
  init_tests.run_tests();

  ViewAccessTests access_tests;
  access_tests.run_tests();

  MiscOperations misc_tests;
  misc_tests.run_tests();

  CopyTests copy_tests;
  copy_tests.run_tests();

  MoveTests move_tests;
  move_tests.run_tests();

  // At this point, we are largely missing tests of subarray_map
  // (which is probably not worth the effor at this moment in time)

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END
