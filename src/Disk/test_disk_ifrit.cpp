/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * SYNOPSIS:
 *
 *    
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
/** 
 *********************************************************************
 *
 * @file      test_hdf5.cpp
 * @brief     Program implementing unit tests for the Hdf5 class
 * @author    James Bordner
 * @date      Thu Feb 21 16:47:35 PST 2008
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include <stdio.h>
#include <string>

#include "hdf5.h"

#include "error.hpp"
#include "test.hpp"
#include "array.hpp"
#include "disk.hpp"

int main(int argc, char ** argv)
{
  Array A;
  A.resize(100,100);

  Scalar * a = A.values();
  int n0,n1;
  A.size(&n0,&n1);

  for (int i0=0; i0<100; i0++) {
    for (int i1=0; i1<100; i1++) {
      int i = i0 + n0*(i1);
      a[i] = i0*3 + i1*5;
    }
  }

  unit_class ("Hdf5");
  unit_open();

  unit_func("file_open");

  Hdf5 hdf5;

  int n[3];

  hdf5.file_open("file_open_test.hdf5","w");
  A.size(&n[0],&n[1],&n[2]);
  hdf5.dataset_open_write ("dataset",n[0],n[1],n[2]);
  hdf5.write(A.values());
  hdf5.dataset_close ();
  hdf5.file_close();

  Array B;

  hdf5.file_open("file_open_test.hdf5","r");
  hdf5.dataset_open_read ("dataset",&n[0],&n[1],&n[2]);
  B.resize(n[0],n[1],n[2]);
  hdf5.read(B.values());
  hdf5.dataset_close ();
  hdf5.file_close();

  Scalar * b = B.values();

  bool passed = true;
  for (int i0=0; i0<100; i0++) {
    for (int i1=0; i1<100; i1++) {
      int i = i0 + n0*(i1);
      if (a[i] != B(i0,i1)) passed = false;
      if (b[i] != A(i0,i1)) passed = false;
    }
  }

  unit_func("read,write");
  unit_assert(passed);

  unit_close();
}
