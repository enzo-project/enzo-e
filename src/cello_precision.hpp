// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef CELLO_PRECISION_HPP
#define CELLO_PRECISION_HPP

/// @file     cello_precision.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Nov 11 17:08:38 PST 2010
/// @brief    Include Cello precision types

enum precision_enum {
  precision_unknown,     //  unknown precision
  precision_default,     //  default precision, based on CONFIG_PRECISION_[SINGLE|DOUBLE]
  precision_single,      //  32-bit field data
  precision_double,      //  64-bit field data
  precision_extended80,  //  80-bit field data
  precision_extended96,  //  96-bit field data
  precision_quadruple,   // 128-bit field data
};

#ifdef CONFIG_PRECISION_SINGLE
#   define default_precision precision_single
typedef float Scalar;
#   define SCALAR_SCANF  "%f"
#   define SCALAR_PRINTF "%e "
#   define SCALAR_MPI     MPI_FLOAT
#   define SCALAR_DEFINED
#endif

#ifdef CONFIG_PRECISION_DOUBLE
#   define default_precision precision_double
typedef double Scalar;
#   define SCALAR_SCANF  "%lf"
#   define SCALAR_PRINTF "%le "
#   define SCALAR_MPI     MPI_DOUBLE
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

  extern const char * precision_name[8];

};

#endif /* CELLO_PRECISION_HPP */
