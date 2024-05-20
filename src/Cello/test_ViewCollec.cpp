// See LICENSE_CELLO file for license and copyright information

/// @file     test_ViewCollec.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2021-11-27
/// @brief    Test program for the ViewCollec

#include "main.hpp"
#include "test.hpp"
#include "view.hpp"

#include "test_ViewTestTools.hpp" // assign_range_3Darray, range_3Darray,
                                  // assert_allequal3D

#include <array>
#include <utility> // for std::pair
#include <vector>

//----------------------------------------------------------------------

namespace {

  template<typename T>
  class ViewCollecFactory {
  public:
    ViewCollecFactory(bool single_array) : single_array_(single_array) { }

    ViewCollec<T> operator()(const std::array<int,3>& shape,
                             const std::vector<std::pair<T,T>>& start_step_list)
      const noexcept
    {
      const std::size_t n_arrays = start_step_list.size();
      if (n_arrays == 0){ return ViewCollec<T>(); }
      ViewCollec<T> collection = init_collec_(n_arrays, shape);

      // now actually initialize the values
      for (std::size_t i = 0; i < n_arrays; i++){
        T start = start_step_list[i].first;
        T step = start_step_list[i].second;
        assign_range_3Darray(collection[i], start, step);
      }
      return collection;
    }

    bool builds_contiguous() const noexcept { return single_array_; }

  private:

    ViewCollec<T> init_collec_(std::size_t n_arrays,
                               const std::array<int,3>& shape) const noexcept
    {
      if (single_array_){
        return ViewCollec<T>(n_arrays, shape);
      } else {
        std::vector<CelloView<T,3>> vec_of_arrays;
        for (std::size_t i = 0; i < n_arrays; i++){
          vec_of_arrays.push_back(CelloView<T,3>(shape[0],shape[1],shape[2]));
        }
        return ViewCollec<T>(vec_of_arrays);
      }
    }

  private:
    bool single_array_;
  };

  template<typename T>
  void common_checks_(const char* test_name, ViewCollec<T>& collec,
                      std::size_t n_arrays, const std::array<int,3>& shape,
                      bool expect_contiguous) noexcept
  {
    ASSERT1(test_name, "The size of the collection is expected to be %zu.",
            n_arrays, collec.size() == n_arrays);

    ASSERT("InitializationTests::test_new_alloc", "Unexpected Shape.",
           (collec.array_shape(0) == shape[0]) &
           (collec.array_shape(1) == shape[1]) &
           (collec.array_shape(2) == shape[2]));

    {
      const char* phrase = (expect_contiguous) ? "should" : "shouldn't";
      ASSERT1(test_name, "The items in the collection %s be contiguous.",
              phrase, collec.contiguous_items() == expect_contiguous);
    }

    if (expect_contiguous){
        CelloView<T,4> backing_arr = collec.get_backing_array();
        ASSERT(test_name, "Unexpected backing_array() shape.",
               (backing_arr.shape(0) == (int)n_arrays) &
               (backing_arr.shape(1) == shape[0]) &
               (backing_arr.shape(2) == shape[1]) &
               (backing_arr.shape(3) == shape[2]));

        bool all_are_aliases = true;
        for (int i = 0; i < static_cast<int>(n_arrays); i++){
          all_are_aliases &= collec[i].is_alias(backing_arr.subarray(i));
        }

        ASSERT(test_name,
               "Each entry of the collection should be a perfect alias of the "
               "corresponding subarray from the underlying backing_array",
               all_are_aliases);
    }
  }

}


//----------------------------------------------------------------------

class InitializationTests{

private:
  void test_default_construction(){
    ViewCollec<double> collec;

    ASSERT("InitializationTests::test_default_construction",
           "a default-constructed collection should have 0 entries.",
           collec.size() == 0);
  }

  // this method tests the copy constructor and move constructor for collec (a
  // pointer to the heap-allocated move-constructed instance is returned by the
  // function).
  //
  // Where applicable it also helps test that the entries of the specified
  // collection alias any specified arrays
  //
  // check_array_elem_values is a functor that accepts a ViewCollec<double>
  // reference and checks 
  // When alias_vec is not empty, this will check that the ith array in collec
  // is an alias of alias_vec
  template<typename F>
  ViewCollec<double>* test_helper_
  (const char* f_name, ViewCollec<double>& collec, F& check_array_elem_values,
   const std::vector<CelloView<double,3>>& alias_vec) const noexcept
  {
    std::size_t n_arrays = collec.size();
    std::array<int,3> shape{collec.array_shape(0),
                            collec.array_shape(1),
                            collec.array_shape(2)};
    bool expect_contiguous = collec.contiguous_items();

    // define a function to verify that the contents of the collection are
    // aliases of the arrays in vec
    auto all_arrays_are_aliases = [=,&alias_vec](ViewCollec<double>& collec)
      {
        if (alias_vec.size() == 0) {
          return true;
        } else if (alias_vec.size() != n_arrays){
          ERROR(f_name,
                "Ill-posed test, alias_vec should either hold no CelloViews "
                "or as many CelloViews as collec");
        }

        for (std::size_t i = 0; i < n_arrays; i++){
          if (!collec[i].is_alias(alias_vec[i])){ return false; }
        }
        return true;
      };

    // first check wrapped aliases for input array
    ASSERT(f_name, ("At least one array from the specified collection "
                    "does not alias an array that it's supposed to"),
           all_arrays_are_aliases(collec));

    // test the copy constructor
    {
      ViewCollec<double> copied_collec(collec);
      common_checks_(f_name, copied_collec, n_arrays, shape, expect_contiguous);
      check_array_elem_values(copied_collec);
      ASSERT(f_name, ("At least one array from the copy-constructed collection "
                      "does not alias an array that it's supposed to"),
             all_arrays_are_aliases(copied_collec));
    }

    // now, test the move constructor (and allocate the destination location)
    ViewCollec<double>* out_ptr = new ViewCollec<double>(std::move(collec));
    common_checks_("InitializationTests::test_wrap_existing", *out_ptr,
                   2, {5,4,3}, expect_contiguous);
    ASSERT(f_name, ("At least one array from the move-constructed collection "
                    "does not alias an array that it's supposed to"),
           all_arrays_are_aliases(*out_ptr));
    check_array_elem_values(*out_ptr);
    return out_ptr;
  }

  /// test construction of a collection that wraps existing arrays
  void test_wrap_existing(){
    ViewCollec<double>* collec_ptr = nullptr; // initialize for later use

    // initialize the reference values
    CelloView<double, 3> ref0 = range_3Darray(5, 4, 3,  1.0, 1.0);
    CelloView<double, 3> ref1 = range_3Darray(5, 4, 3, -1.0,-1.0);
    auto check_array_elem_values = [&ref0, &ref1](ViewCollec<double>& collec)
      {
        assert_allequal3D(collec[0], ref0);
        assert_allequal3D(collec[1], ref1);
      };


    {
      // define the vector of arrays that will be wrapped
      std::vector<CelloView<double, 3>> vec{ref0.deepcopy(), ref1.deepcopy()};

      // initialize the collection
      ViewCollec<double> temp_collec(vec);

      // test properties of temp_collec
      common_checks_("InitializationTests::test_wrap_existing", temp_collec,
                     2, {5,4,3}, false);
      check_array_elem_values(temp_collec);

      // check that temp_collec actually wraps the entries in vec, and then
      // test the copy & move constructors
      collec_ptr = test_helper_("InitializationTests::test_wrap_existing",
                                temp_collec, check_array_elem_values, vec);
      // the contents of temp_collec have been moved to collec_ptr
    }

    // finally, repeat some of the checks now that the wrapped vector of arrays
    // and the variable that originally held the collection are out of scope
    // (to confirm there are no ref-counting issues)
    common_checks_("InitializationTests::test_wrap_existing", *collec_ptr,
                   2, {5,4,3}, false);
    check_array_elem_values(*collec_ptr);

    delete collec_ptr;
  }

  /// test construction of a collection that allocates it's own arrays
  void test_new_alloc(){
    ViewCollec<double>* collec_ptr = nullptr; // initialize for later use

    // initialize the reference values
    CelloView<double, 3> ref0 = range_3Darray(5, 4, 3,  1.0, 1.0);
    CelloView<double, 3> ref1 = range_3Darray(5, 4, 3, -1.0,-1.0);
    auto check_array_elem_values = [&ref0, &ref1](ViewCollec<double>& collec)
      {
        assert_allequal3D(collec[0], ref0);
        assert_allequal3D(collec[1], ref1);
      };

    {
      // initialize the collection
      ViewCollec<double> temp_collec(2, {5,4,3});

      // test properties of temp_collec
      common_checks_("InitializationTests::test_new_alloc", temp_collec,
                     2, {5,4,3}, true);

      // update the values within the array
      ref0.copy_to(temp_collec[0]);
      ref1.copy_to(temp_collec[1]);

      // test properties of temp_collec
      common_checks_("InitializationTests::test_new_alloc", temp_collec,
                     2, {5,4,3}, true);
      check_array_elem_values(temp_collec);

      CkPrintf("Testing the copy/move constructors\n");
      fflush(stdout);

      // test the copy & move constructors
      collec_ptr = test_helper_("InitializationTests::test_new_alloc",
                                temp_collec, check_array_elem_values, {});
      // the contents of temp_collec have been moved to collec_ptr
    }

    // finally, repeat some of the checks now that the vector of arrays left
    // the scope (to confirm there are no ref-counting issues)
    common_checks_("InitializationTests::test_new_alloc", *collec_ptr,
                   2, {5,4,3}, true);
    check_array_elem_values(*collec_ptr);

    delete collec_ptr;
  }

  /// this is more of a test of the factory methods than anything else
  void test_with_factories(){
    for (bool use_single_arr: {false, true}){
      ViewCollecFactory<double> factory(use_single_arr);
      ViewCollec<double> collec = factory({3,4,5},
                                          {{1.,0.5}, {100., 2.}, {1000., 4.}});

      common_checks_("test_with_factories", collec, 3, {3,4,5}, use_single_arr);
    }
  }

public:

  void run_tests(){
    test_default_construction();
    test_wrap_existing();
    test_new_alloc();
    test_with_factories();
  }

};

//----------------------------------------------------------------------

class VariableAssignmentTests{
  // these are tests that check that check the assignments of arrays to
  // variables.

  template<bool test_copy_assignment>
  void test_assignment_helper_
  (const ViewCollecFactory<double>& main_factory,
   const ViewCollecFactory<double>& alt_factory,
   const std::array<int,3>& alt_shape,
   const std::vector<std::pair<double,double>>& alt_start_step_list)
  {
    std::string test_name;

    // initialize the destination variable
    ViewCollec<double> dest_collec(alt_factory(alt_shape, alt_start_step_list));

    // perform the assignment
    if (test_copy_assignment){
      test_name = "VariableAssignmentTests::test_assignment_helper_<true>";
      // make the source variable a constant so that the compiler cannot use
      // move assignment
      const ViewCollec<double> src(main_factory({3,4,5},
                                                {{1.,0.5}, {4.,2.}, {16.,4.}}));
      dest_collec = src;
    } else {
      test_name = "VariableAssignmentTests::test_assignment_helper_<false>";
      ViewCollec<double> src(main_factory({3,4,5},
                                          {{1.,0.5}, {4.,2.}, {16.,4.}}));
      dest_collec = std::move(src);
    }

    bool expect_contiguous = main_factory.builds_contiguous();
    common_checks_(test_name.c_str(), dest_collec, 3, {3,4,5},
                   expect_contiguous);
  }

  void test_assignment_(const ViewCollecFactory<double>& main_factory){
    ViewCollecFactory<double> contig_factory(true);
    ViewCollecFactory<double> wrapped_arr_factory(false);

    // test copy assignment
    test_assignment_helper_<true>(main_factory, contig_factory,
                                  {1,2,3}, {{-1., -0.5}, {-4., -2.}});
    test_assignment_helper_<true>(main_factory, wrapped_arr_factory,
                                  {3,2,1}, {{-1., -0.5}, {-4., -2.}});

    // test move assignment
    test_assignment_helper_<false>(main_factory, contig_factory,
                                   {1,2,3}, {{-1., -0.5}, {-4., -2.}});
    test_assignment_helper_<false>(main_factory, wrapped_arr_factory,
                                   {3,2,1}, {{-1., -0.5}, {-4., -2.}});
  }

public:

  void run_tests(){
    test_assignment_(ViewCollecFactory<double>(true));
    test_assignment_(ViewCollecFactory<double>(false));
  }

};

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{
  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ViewCollec");
  
  InitializationTests init_tests;
  init_tests.run_tests();

  VariableAssignmentTests var_assign_tests;
  var_assign_tests.run_tests();

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END
