//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Error and warning definitions

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

#define WARNING(MESSAGE) { \
   printf ("WARNING %s:%d %s\n",__FILE__,__LINE__,MESSAGE);	\
  fflush(stdout); \
}

#define ERROR(MESSAGE) { \
   printf ("ERROR %s:%d %s\n",__FILE__,__LINE__,MESSAGE); \
  fflush(stdout); \
}

#define NOT_IMPLEMENTED(X) { \
   printf ("%s:%d WARNING: " X " is not implemented yet\n",__FILE__,__LINE__); \
   fflush(stdout); \
}
