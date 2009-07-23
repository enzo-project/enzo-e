//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      test_parameters.cpp
 * @brief     Program implementing unit tests for the Parameters class
 * @author    James Bordner
 * @date      Thu Feb 21 16:04:03 PST 2008
 * 
 * $Id: test_parameters.cpp 715 2009-07-08 23:48:09Z bordner $
 * 
 *********************************************************************
 */
 
#include <stdio.h>

#include "test.hpp"
#include "error.hpp"
#include "parameters.hpp"

main()
{

  unit_class ("Parameters");

  unit_open();

  //----------------------------------------------------------------------
  // test parameter
  //----------------------------------------------------------------------

  try {
    Parameters parameters;

    FILE * file_pointer = fopen ("in.test_parameters","r");
    parameters.read(file_pointer);

    bool did_pass = (parameters.get_value("key1") == "value1");
    
    unit_func("get_value()");
    unit_assert(did_pass);

  }
  // ERRORS ARE NOT CAUGHT HERE--WHY?
  catch (ExceptionBadPointer) {
    printf ("CELLO ERROR: Bad pointer.\n");
  }
  catch (...) {
    printf ("CELLO ERROR: Unknown error.\n",error_message);
  }

  unit_close();

}
