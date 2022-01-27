
#ifdef FORTRAN
#  ifdef CONFIG_PRECISION_SINGLE
#     define CONFIG_BFLOAT_4
#     define tiny 1.e-20
#     define huge 1.e+20
#     define P_PREC real (kind=4)
#     define R_PREC real (kind=4)
#     define CMPLX_PREC complex(kind=8)
      integer, parameter :: RKIND=4
      integer, parameter :: PKIND=4
#  endif
#   ifdef CONFIG_PRECISION_DOUBLE
#     define CONFIG_BFLOAT_8
#     define tiny 1.d-20
#     define huge 1.d+20
#     define P_PREC real(kind=8)
#     define R_PREC real(kind=8)
#     define CMPLX_PREC complex(kind=16)
      integer, parameter :: PKIND=8
      integer, parameter :: RKIND=8
#   endif
#   ifdef CONFIG_PRECISION_QUAD
#     define tiny 1.d-35
#     define huge 1.d+35
#     define P_PREC real(kind=8)
      integer, parameter :: PKIND=8
      integer, parameter :: RKIND=8
#   endif
#   ifdef SMALL_INTS
#     define INTG_PREC integer (kind=4)
#     define LOGIC_PREC logical (kind=4)
      integer, parameter :: IKIND=4
      integer, parameter :: LKIND=4
#   endif
#   ifdef LARGE_INTS
#     define INTG_PREC integer (kind=8)
#     define LOGIC_PREC logical (kind=8)
      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8
#   endif
#endif
#define ERROR_CALC_DT_NAN                 100
#define ERROR_CIC_DEPOSIT_CLOUD_SIZE      200
#define ERROR_EULER_GE_LT_0_A             300
#define ERROR_EULER_GE_LT_0_B             400
#define ERROR_EULER_E_LT_0                500
#define ERROR_FLUX_HLL_E_LT_0_A           600
#define ERROR_FLUX_HLL_E_LT_0_B           700
#define ERROR_FLUX_HLL_E_LT_0_C           800
#define ERROR_FLUX_HLL_E_LT_0_D           900
#define ERROR_FLUX_HLL_COL_LT_0          1000
#define ERROR_FLUX_TWOSHOCK_E_LT_0_A     1100
#define ERROR_FLUX_TWOSHOCK_E_LT_0_B     1200
#define ERROR_FLUX_TWOSHOCK_E_LT_0_C     1300
#define ERROR_FLUX_TWOSHOCK_E_LT_0_D     1400
#define ERROR_FLUX_TWOSHOCK_COL_LT_0     1500
#define ERROR_INTERP1D_REFINE_MAX        1600
#define ERROR_INTERP2D_REFINE_MAX        1700
#define ERROR_INTERP3D_REFINE_MAX        1800
#define ERROR_INTERP3D_REFINE_PARENT_NAN 1900
#define ERROR_INTERP3D_REFINE_FRAC       2000
#define ERROR_INTERP3D_REFINE_NNEG       2100
#define ERROR_INTERP3D_REFINE_S_NAN      2200
#define ERROR_INTERPOLATE_NDIM           2300
#define ERROR_INTERPOLATE_METHOD_RANGE   2400
#define ERROR_INTERPOLATE_BOUNDS         2500
#define ERROR_INTERPOLATE_INDICES        2600
#define ERROR_INTERPOLATE_METHOD_2       2700
#define ERROR_LIN3D_REFINE               2800
#define ERROR_LIN3D_NEG                  2900
#define ERROR_PGAS2D_DUAL_GE_LT_0        3000
#define ERROR_PPM_DE_DIM                 3100

#define ENZO_ERROR_EULER_EU1                   10101
#define ENZO_ERROR_EULER_GESLICE_1             10102
#define ENZO_ERROR_EULER_GESLICE_2             10103
#define ENZO_ERROR_FLUX_HLL_COLOR              10201
#define ENZO_ERROR_FLUX_HLL_DSLICE             10202
#define ENZO_ERROR_FLUX_HLL_ESLICE             10203
#define ENZO_ERROR_FLUX_TWOSHOCK_1             10301
#define ENZO_ERROR_FLUX_TWOSHOCK_DSLICE        10302
#define ENZO_ERROR_FLUX_TWOSHOCK_ESLICE        10303
#define ENZO_ERROR_FLUX_TWOSHOCK_FALLBACK      10304
#define ENZO_ERROR_INTERPOLATE_METHOD          10401
#define ENZO_ERROR_INTERPOLATE_NDIM            10402
#define ENZO_ERROR_INTERPOLATE_PSTART_IMETHOD2 10403
#define ENZO_ERROR_INTERPOLATE_PSTART_PEND     10404
#define ENZO_ERROR_INTERPOLATE_PSTART_REFINE   10405
#define ENZO_ERROR_XEULER_DSLICE               10501
#define ENZO_ERROR_XEULER_ESLICE               10502
#define ENZO_ERROR_XEULER_GESLICE              10503
#define ENZO_ERROR_YEULER_DSLICE               10601
#define ENZO_ERROR_YEULER_ESLICE               10602
#define ENZO_ERROR_YEULER_GESLICE              10603
#define ENZO_ERROR_ZEULER_DSLICE               10701
#define ENZO_ERROR_ZEULER_ESLICE               10702
#define ENZO_ERROR_ZEULER_GESLICE              10703
