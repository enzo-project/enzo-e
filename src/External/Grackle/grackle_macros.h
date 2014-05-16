/***********************************************************************
/
/ Grackle definitions
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_MACROS_H_
#define __GRACKLE_MACROS_H_
/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

#ifdef CONFIG_THROW_ABORT
#define GRACKLE_FAIL(A) raise(SIGABRT);
#define GRACKLE_VFAIL(A, ...) raise(SIGABRT);
#else
#define GRACKLE_FAIL(A) throw(GrackleFatalException(A, __FILE__, __LINE__));
#define GRACKLE_VFAIL(format, ...) {snprintf(current_error, 254, format, ##__VA_ARGS__); throw(GrackleFatalException(current_error, __FILE__, __LINE__));}
#endif

#define FORTRAN_NAME(NAME) NAME##_

/* HDF5 definitions */

#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I8 H5T_STD_I64BE
#define HDF5_FILE_R4 H5T_IEEE_F32BE
#define HDF5_FILE_R8 H5T_IEEE_F64BE
#define HDF5_FILE_B8 H5T_STD_B8BE

#define HDF5_I4  H5T_NATIVE_INT
#define HDF5_I8  H5T_NATIVE_LLONG
#define HDF5_R4  H5T_NATIVE_FLOAT
#define HDF5_R8  H5T_NATIVE_DOUBLE
#define HDF5_R16 H5T_NATIVE_LDOUBLE

/* Precision-dependent definitions */

#define LARGE_INTS

#ifdef CONFIG_PRECISION_SINGLE
#  define CONFIG_BFLOAT_4
#endif
#ifdef CONFIG_PRECISION_DOUBLE
#  define CONFIG_BFLOAT_8
#endif


#ifdef SMALL_INTS
#define ISYM "d"
#define nint(A) ( (int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long long) ((A) + 0.5*sign(A)) )
#define ABS(A) abs((int) (A))
#endif

#ifdef LARGE_INTS
#define ISYM "lld"
#define nint(A) ( (long long) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long long) ((A) + 0.5*sign(A)) )
#define ABS(A) labs((long long) (A))
#endif


#ifdef CONFIG_BFLOAT_4
#define FSYM "f"
#define ESYM "e"
#endif

#ifdef CONFIG_BFLOAT_8
#define FSYM "lf"
#define ESYM "le"
#endif

#define GSYM "g"

/* Standard definitions (well, fairly standard) */

#ifndef NULL
#define NULL      0
#endif

#ifdef FAIL
#undef FAIL
#endif
#define FAIL      0
#define SUCCESS   1

#ifndef FALSE
#define FALSE     0
#define TRUE      1
#endif

#define FLOAT_UNDEFINED  -99999.0
#define INT_UNDEFINED    -99999
#define MAX_LINE_LENGTH                   512

/* Macro definitions (things C should have) */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define POW(X,Y) pow((double) (X), (double) (Y))
#define COS(X) cos((double) (X))
#define SIN(X) sin((double) (X))

#endif
