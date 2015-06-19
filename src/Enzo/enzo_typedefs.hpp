#ifndef ENZO_TYPEDEFS_HPP
#define ENZO_TYPEDEFS_HPP

//----------------------------------------------------------------------
// TYPEDEFS
//----------------------------------------------------------------------

typedef int                 Eint32;     // c_message only
typedef long long           long_int;   // use long long
typedef long long int       Elong_int;  // use long long
typedef long long unsigned  global_index; // 

typedef long                gr_int;  // Grackle int

#ifdef CONFIG_PRECISION_SINGLE
   typedef float       enzo_float;
   typedef float       gr_float;     // Grackle float
#elif  CONFIG_PRECISION_DOUBLE
   typedef double      enzo_float;
   typedef double      gr_float;     // Grackle float
#elif  CONFIG_PRECISION_QUAD
   typedef long double enzo_float;
   typedef long double gr_float;     // Grackle float
#else
#  error "Must define CONFIG_PRECISION_[SINGLE|DOUBLE|QUAD]"
#endif

/* #include "enzo_typedefs_30.hpp" */

#endif /* ENZO_TYPEDEFS_HPP */

