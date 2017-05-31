#ifdef FORTRAN
#ifdef CONFIG_PRECISION_SINGLE
      integer, parameter :: RKIND=4
#endif

#ifdef CONFIG_PRECISION_DOUBLE
      integer, parameter :: RKIND=8
#endif

#ifdef CONFIG_PRECISION_QUAD
      integer, parameter :: RKIND=8
#endif
#endif
