//345678901234567890123456789012345678901234567890123456789012345678901234567890

#define _TRACE_ if (trace) { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); printf ("TRACE %d %s:%d\n",_ip_,__FILE__,__LINE__); fflush(stdout); }
#define _BARRIER_ if (trace) { MPI_Barrier (MPI_COMM_WORLD); }

#define _TEMPORARY_ { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); printf ("TEMPORARY %d %s:%d\n",_ip_,__FILE__,__LINE__); fflush(stdout); }
#define _WARNING_(X) { int _ip_; MPI_Comm_rank (MPI_COMM_WORLD, &_ip_); if (_ip_==0) printf ("WARNING %d %s:%d" X "\n",__FILE__,__LINE__); fflush(stdout); }

#define NOT_IMPLEMENTED(X) printf ("%s:%d WARNING: " X " is not implemented yet\n",__FILE__,__LINE__);

#define MIN(a,b) ( (a) < (b) ? (a) : (b))
#define MAX(a,b) ( (a) > (b) ? (a) : (b))



