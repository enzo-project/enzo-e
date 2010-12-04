// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef CELLO_PRECISION_HPP
#define CELLO_PRECISION_HPP

/// @file     cello_precisionhpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Nov 11 17:08:38 PST 2010
/// @brief    Include Cello precision types

enum precision_enum {
  precision_unknown,     //  unknown precision
  precision_default,     //  default precision, based on CONFIG_PRECISION_[SINGLE|DOUBLE]
  precision_half,       //   16-bit field data
  precision_single,      //  32-bit field data
  precision_double,      //  64-bit field data
  precision_extended80,  //  80-bit field data
  precision_extended96,  //  96-bit field data
  precision_quadruple,   // 128-bit field data
};

#ifdef CONFIG_PRECISION_SINGLE
#   define default_precision precision_single
#   define Scalar float
#   define SCALAR_SCANF  "%f"
#   define SCALAR_PRINTF "%e "
#   define SCALAR_MPI     MPI_FLOAT
#   define SCALAR_DEFINED
#endif

#ifdef CONFIG_PRECISION_DOUBLE
#   define default_precision precision_double
#   define Scalar double
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

  int precision_size     (enum precision_enum);
  int precision_supported(enum precision_enum);

  extern const char * precision_name[8];

};

#endif /* CELLO_PRECISION_HPP */
