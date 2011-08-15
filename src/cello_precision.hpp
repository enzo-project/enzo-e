// See LICENSE_CELLO file for license and copyright information

#ifndef CELLO_PRECISION_HPP
#define CELLO_PRECISION_HPP

/// @file     cello_precision.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Nov 11 17:08:38 PST 2010
/// @brief    Include Cello precision types

enum precision_enum {
  // @@@ KEEP IN SYNCH WITH precision_name in cello_precision.cpp
  precision_unknown,     //  unknown precision
  precision_default,     //  default precision
  precision_single,      //  32-bit field data
  precision_double,      //  64-bit field data
  precision_extended80,  //  80-bit field data
  precision_extended96,  //  96-bit field data
  precision_quadruple,   // 128-bit field data
};

#ifdef CONFIG_PRECISION_SINGLE
#   define default_precision precision_single
#   define SCALAR_DEFINED
#endif

#ifdef CONFIG_PRECISION_DOUBLE
#   define default_precision precision_double
#   ifdef SCALAR_DEFINED
#      error Both CONFIG_PRECISION_SINGLE and CONFIG_PRECISION_DOUBLE defined
#   endif
#   define SCALAR_DEFINED
#endif

#ifndef SCALAR_DEFINED
#   error Neither CONFIG_PRECISION_SINGLE nor CONFIG_PRECISION_DOUBLE defined
#endif

namespace cello {

  int sizeof_precision       (enum precision_enum);
  int is_precision_supported (enum precision_enum);
  double machine_epsilon     (enum precision_enum);
  extern const char * precision_name[7];

  template <class T>
  T err_rel (T a, T b)
  {  return (a != 0.0) ? fabs((a - b) / a) : fabs(a-b);  }

  template <class T>
  T err_abs (T a, T b)
  {  return fabs(a-b);  };

};

#endif /* CELLO_PRECISION_HPP */
