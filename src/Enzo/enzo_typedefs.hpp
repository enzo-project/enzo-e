#ifndef ENZO_TYPEDEFS_HPP
#define ENZO_TYPEDEFS_HPP

//----------------------------------------------------------------------
// TYPEDEFS
//----------------------------------------------------------------------

typedef int                 Eint32;     // c_message only
typedef long long           long_int;   // use long long
typedef long long int       Elong_int;  // use long long
typedef long long unsigned  global_index; // 

#ifdef CONFIG_PRECISION_SINGLE
   typedef float       enzo_float;
#elif  CONFIG_PRECISION_DOUBLE
   typedef double      enzo_float;
#elif  CONFIG_PRECISION_QUAD
   typedef long double enzo_float;
#else
#  error "Must define one of CONFIG_PRECISION_[SINGLE|DOUBLE|QUAD]"
#endif

// we define these so that we can reduce the number of ifdefs scattered
// throughout the codebase (in the future, it probably makes sense to remove
// these datatypes from all areas of the codebase other than the chemistry
// component).
#ifndef CONFIG_USE_GRACKLE
   typedef long       gr_int;
   typedef enzo_float gr_float;
#endif

/// alias for EFlt3DArray
typedef CelloView<enzo_float,3> EFlt3DArray;

typedef std::vector<std::string> str_vec_t;

/* #include "enzo_typedefs_30.hpp" */

#endif /* ENZO_TYPEDEFS_HPP */

