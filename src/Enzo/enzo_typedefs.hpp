#ifndef ENZO_TYPEDEFS_H
#define ENZO_TYPEDEFS_H

//----------------------------------------------------------------------
// TYPEDEFS
//----------------------------------------------------------------------

typedef int            Eint32;     // c_message only
typedef long long      long_int;   // use long long
typedef long long int  Elong_int;  // use long long
typedef long long unsigned  global_index; // 

#ifdef CONFIG_PRECISION_SINGLE
   typedef float enzo_float;
#elif  CONFIG_PRECISION_DOUBLE
   typedef double enzo_float;
#elif  CONFIG_PRECISION_QUADRUPLE
   typedef long double enzo_float;
#else
#  error "Must define CONFIG_PRECISION_[SINGLE|DOUBLE|QUADRUPLE]"
#endif


#endif /* ENZO_TYPEDEFS_H */

