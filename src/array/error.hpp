//345678901234567890123456789012345678901234567890123456789012345678901234567890

#ifndef ERROR_HPP
#define ERROR_HPP

// $Id$
/**
 * @file    error.hpp
 * @brief   Error and warning definitions
 * @author  James Bordner 
 * @version 1.0
 *
 * Definitions
 *
 * (*) WARNING(MESSAGE)
 * (*) ERROR(MESSAGE)
 * (*) NOT_IMPLEMENTED(X)
 *
 */
// $Log$

/// Handle a (non-lethal) warning message

#define WARNING(MESSAGE) { \
   printf ("WARNING %s:%d %s\n",__FILE__,__LINE__,MESSAGE);	\
  fflush(stdout); \
}

/// Handle a (lethal) error message

#define ERROR(MESSAGE) { \
   printf ("ERROR %s:%d %s\n",__FILE__,__LINE__,MESSAGE); \
  fflush(stdout); \
  exit(1); // Need more graceful shutdown, maybe via Control::exit() \
}

/// Placeholder for functions that are not implemented yet

#define NOT_IMPLEMENTED(X) { \
   printf ("%s:%d WARNING: " X " is not implemented yet\n",__FILE__,__LINE__); \
   fflush(stdout); \
}

#endif
