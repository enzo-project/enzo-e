#define _TRACE_ if (trace) { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); printf ("TRACE %d %s:%d\n",_ip_,__FILE__,__LINE__); fflush(stdout); }

#define _TEMPORARY_ { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); printf ("TEMPORARY %d %s:%d\n",_ip_,__FILE__,__LINE__); fflush(stdout); }

#define NOT_IMPLEMENTED(X) printf ("%s:%d WARNING: " X " is not implemented yet\n",__FILE__,__LINE__);

#define MIN(a,b) ( (a) < (b) ? (a) : (b))
#define MAX(a,b) ( (a) > (b) ? (a) : (b))



