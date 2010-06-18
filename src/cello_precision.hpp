
enum precision_type {
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
#   define SCALAR_HDF5    H5T_NATIVE_FLOAT
#   define SCALAR_DEFINED
#endif

#ifdef CONFIG_PRECISION_DOUBLE
#   define default_precision precision_double
#   define Scalar double
#   define SCALAR_SCANF  "%lf"
#   define SCALAR_PRINTF "%le "
#   define SCALAR_MPI     MPI_DOUBLE
#   define SCALAR_HDF5    H5T_NATIVE_DOUBLE
#   define SCALAR_DEFINED /* if error here, then both single and double defined */
#endif

#ifndef SCALAR_DEFINED
#   error CONFIG_PRECISION_* not defined
#endif

namespace cello {

  int precision_size     (enum precision_type);
  int precision_supported(enum precision_type);
  int precision_supported(enum precision_type);
  long double precision_array_value 
  ( void * array, int index, enum precision_type precision );

  extern const char * precision_name[8];

};


