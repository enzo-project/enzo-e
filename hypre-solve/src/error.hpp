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

#define WARNING(MESSAGE) { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); printf ("TRACE %d %s:%d\n",_ip_,__FILE__,__LINE__); fflush(stdout); }
