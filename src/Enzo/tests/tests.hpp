// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/tests/tests.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for tests subcomponent within the \ref Enzo component
///
/// This subcomponent contains initializers for setting up test-problems

#ifndef ENZO_TESTS_TESTS_HPP
#define ENZO_TESTS_TESTS_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/mesh.hpp"    // Block
#include "Cello/problem.hpp" // Initial

#include "Enzo/enzo.hpp"

//----------------------------------------------------------------------
// Component headers
//----------------------------------------------------------------------

#include "tests/EnzoInitialFeedbackTest.hpp"
#include "tests/EnzoInitialGrackleTest.hpp"
#include "tests/EnzoInitialPpmlTest.hpp"
#include "tests/EnzoInitialMergeSinksTest.hpp"
#include "tests/EnzoInitialAccretionTest.hpp"
#include "tests/EnzoInitialBBTest.hpp"

#endif /* ENZO_TESTS_TESTS_HPP */
