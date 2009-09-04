//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2009 James Bordner
 * Copyright (C) 2009 Laboratory for Computational Astrophysics
 * Copyright (C) 2009 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      test_memory.cpp
 * @brief     Program implementing unit tests for the Memory component
 * @author    James Bordner
 * @date      Thu Sep  3 16:37:24 PDT 2009
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include <stdio.h>
#include <string>

#include "error.hpp"

#include "memory.hpp"

#include "test.hpp"

main()
{
  unit_class ("Memory");
  unit_open();

  unit_func("allocate");
  double * a_double = (double *) memory.allocate(sizeof(double[10]));
  unit_assert(false);

  unit_close();
}
