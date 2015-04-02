// See LICENSE_CELLO file for license and copyright information

/// @file     gtest_Colormap.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the ColormapRGB class

#include "main.hpp"
#include "test.hpp"

#include "io.hpp"

#include "gtest/gtest.h"

class Foo { public: Foo() {};
  int Bar (std::string a, std::string b)
  {return 0;}
};

namespace {

// The fixture for testing class Foo.
class FooTest : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  FooTest() {
    // You can do set-up work for each test here.
  }

  virtual ~FooTest() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Objects declared here can be used by all tests in the test case for Foo.
};

// Tests that the Foo::Bar() method does Abc.
TEST_F(FooTest, MethodBarDoesAbc) {
  const std::string input_filepath = "this/package/testdata/myinputfile.dat";
  const std::string output_filepath = "this/package/testdata/myoutputfile.dat";
  Foo f;
  EXPECT_EQ(0, f.Bar(input_filepath, output_filepath));
}

// Tests that Foo does Xyz.
TEST_F(FooTest, DoesXyz) {
  // Exercises the Xyz feature of Foo.
}

}  // namespace

PARALLEL_MAIN_BEGIN
{
  ::testing::InitGoogleTest(&PARALLEL_ARGC, PARALLEL_ARGV);
  RUN_ALL_TESTS();
  PARALLEL_EXIT;
}
PARALLEL_MAIN_END
