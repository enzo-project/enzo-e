
#ifdef FORTRAN
#  ifdef CONFIG_PRECISION_SINGLE
#     define tiny 1.d-20
#     define huge 1.d+20
#     define P_PREC real*4
#     define R_PREC real*4
#     define CMPLX_PREC complex*8
      integer, parameter :: RKIND=4
      integer, parameter :: PKIND=4
#  endif

#   ifdef CONFIG_PRECISION_DOUBLE
#     define tiny 1.d-20
#     define huge 1.d+20
#     define P_PREC real*8
#     define R_PREC real*8
#     define CMPLX_PREC complex*16
      integer, parameter :: PKIND=8
      integer, parameter :: RKIND=8
#   endif

#   ifdef CONFIG_PRECISION_QUAD
#     define tiny 1.d-35
#     define huge 1.d+35
#     define P_PREC real*8
      integer, parameter :: PKIND=8
      integer, parameter :: RKIND=8
#   endif

#   ifdef SMALL_INTS
#     define INTG_PREC integer*4
#     define LOGIC_PREC logical*4
      integer, parameter :: IKIND=4
      integer, parameter :: LKIND=4
#   endif
	
#   ifdef LARGE_INTS
#     define INTG_PREC integer*8
#     define LOGIC_PREC logical*8
      integer, parameter :: IKIND=8
      integer, parameter :: LKIND=8
#   endif
	
#endif
