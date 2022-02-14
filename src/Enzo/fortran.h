#include "enzo_defines.hpp"

#ifdef CONFIG_PRECISION_SINGLE
#  define ENZO_REAL real (kind=4)
#  define ENZO_TINY 1e-10
#endif

#ifdef CONFIG_PRECISION_DOUBLE
#  define ENZO_REAL real (kind=8)
#  define ENZO_TINY 1d-20
#endif

#ifdef CONFIG_PRECISION_QUAD
#  define ENZO_REAL real (kind=16)
#  define ENZO_TINY 1d-40
#endif

