//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Performance type definitions

/**
 * 
 * @file      performance.hpp
 * @brief     lcaperf-related defines
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 */

#ifdef USE_LCAPERF

#   include "lcaperf.h"


#   define LCAPERF_BEGIN(SEGMENT) lcaperf.begin(SEGMENT)
#   define LCAPERF_END(SEGMENT)   lcaperf.end(SEGMENT)
#   define LCAPERF_START(REGION)  lcaperf.start(REGION)
#   define LCAPERF_STOP(REGION)   lcaperf.stop(REGION)
#   define LCAPERF_GLOBAL(NAME,VALUE) lcaperf.global(NAME,VALUE)

#else

#   define LCAPERF_BEGIN(SEGMENT)     ;
#   define LCAPERF_END(SEGMENT)       ;
#   define LCAPERF_START(REGION)      ;
#   define LCAPERF_STOP(REGION)       ;
#   define LCAPERF_GLOBAL(NAME,VALUE) ;

#endif



