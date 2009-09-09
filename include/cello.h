#ifndef CELLO_DEF
#define CELLO_DEF

/*********************************************************************
* Define either CONFIG_PRECISION_SINGLE or CONFIG_PRECISION_DOUBLE
**********************************************************************/

#define CONFIG_PRECISION_SINGLE
/* #define CONFIG_PRECISION_DOUBLE */

/*********************************************************************
* GLOBAL DECLARATIONS
**********************************************************************/

#ifdef CONFIG_PRECISION_SINGLE
#   define Scalar float
#   define SCALAR_SCANF  "%f"
#   define SCALAR_PRINTF "%e "
#   define SCALAR_MPI     MPI_FLOAT
#   define SCALAR_HDF5    H5T_NATIVE_FLOAT
#endif

#ifdef CONFIG_PRECISION_DOUBLE
#   define Scalar double
#   define SCALAR_SCANF  "%lf"
#   define SCALAR_PRINTF "%le "
#   define SCALAR_MPI    MPI_DOUBLE
#   define SCALAR_HDF5    H5T_NATIVE_DOUBLE
#endif

#endif /* SCALAR_DEF */

#endif /* CELLO_DEF */
