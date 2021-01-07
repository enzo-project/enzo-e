#include "enzo_defines.hpp"

#define DEBUG_MESSAGE CALL f_warning(__FILE__,__LINE__)
#define ERROR_MESSAGE CALL f_error(__FILE__,__LINE__)
#define WARNING_MESSAGE CALL f_warning(__FILE__,__LINE__)

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

