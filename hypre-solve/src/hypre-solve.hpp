#define _TRACE_ if (trace) { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); printf ("%d %s:%d\n",_ip_,__FILE__,__LINE__); fflush(stdout); }

#define MIN(a,b) ( (a) < (b) ? (a) : (b))
#define MAX(a,b) ( (a) > (b) ? (a) : (b))

