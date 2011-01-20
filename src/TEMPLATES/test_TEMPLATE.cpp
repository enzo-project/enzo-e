// $Id: test_Template.cpp 1696 2010-08-04 05:56:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_TEMPLATE.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Template class

#include "cello.hpp"


#include "error.hpp"
#include "test.hpp"
#include "class.hpp"
#include "field.hpp"

#include "parallel.def"
#include PARALLEL_CHARM_INCLUDE(test_Template.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  unit_class ("TemplateDescr");

  Template template;

  unit_assert (false)

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Template.def.h)
