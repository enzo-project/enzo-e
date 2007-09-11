
/// Performance type definitions

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

#ifdef USE_JBPERF

#   include "jbPerf.h"


#   define JBPERF_BEGIN(SEGMENT) jbPerf.begin(SEGMENT)
#   define JBPERF_END(SEGMENT)   jbPerf.end(SEGMENT)
#   define JBPERF_START(REGION)  jbPerf.start(REGION)
#   define JBPERF_STOP(REGION)   jbPerf.stop(REGION)
#   define JBPERF_GLOBAL(NAME,VALUE) jbPerf.global(NAME,VALUE)

#else

#   define JBPERF_BEGIN(SEGMENT)     ;
#   define JBPERF_END(SEGMENT)       ;
#   define JBPERF_START(REGION)      ;
#   define JBPERF_STOP(REGION)       ;
#   define JBPERF_GLOBAL(NAME,VALUE) ;

#endif



