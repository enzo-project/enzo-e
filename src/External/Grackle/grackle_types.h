/***********************************************************************
/
/ Grackle variable types
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_TYPES_H_
#define __GRACKLE_TYPES_H_
/***********************************************************************
/  
/ VARIABLE TYPES
/
************************************************************************/

#ifdef SMALL_INTS
#define gr_int int
#endif

#ifdef LARGE_INTS
#define gr_int long long
#endif

#ifdef CONFIG_BFLOAT_4
#define gr_float float
#endif

#ifdef CONFIG_BFLOAT_8
#define gr_float double
#endif

#endif
