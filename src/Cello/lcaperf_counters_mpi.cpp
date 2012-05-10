// $Id: counters_mpi.cpp 2093 2011-03-12 01:17:05Z bordner $
// See LICENSE file for license and copyright information

/// @file     counters_mpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-22
/// @brief    Implementation of the CountersMpi class

#ifdef CONFIG_USE_MPI
#include "performance.hpp"
namespace lca {

double CountersMpi::mpi_time_start_ = 0;
long long CountersMpi::mpi_call_time_ = 0;

long long CountersMpi::mpi_calls_      = 0;
long long CountersMpi::mpi_time_       = 0;

long long CountersMpi::mpi_send_time_ = 0;
long long CountersMpi::mpi_send_bytes_ = 0;
long long CountersMpi::mpi_send_calls_ = 0;

long long CountersMpi::mpi_recv_time_ = 0;
long long CountersMpi::mpi_recv_bytes_ = 0;
long long CountersMpi::mpi_recv_calls_ = 0;

long long CountersMpi::mpi_sync_time_ = 0;
long long CountersMpi::mpi_sync_procs_ = 0;
long long CountersMpi::mpi_sync_calls_ = 0;

//----------------------------------------------------------------------

CountersMpi::CountersMpi() throw ()
  : CountersUser(),
    ip_mpi_(-1),
    np_mpi_(-1),
    ip_node_(-1),
    np_node_(-1)
{
  create("mpi-time",       counter_type_relative);
  create("mpi-calls",      counter_type_relative);
  create("mpi-send-time",  counter_type_relative);
  create("mpi-send-bytes", counter_type_relative);
  create("mpi-send-calls", counter_type_relative);
  create("mpi-recv-time",  counter_type_relative);
  create("mpi-recv-bytes", counter_type_relative);
  create("mpi-recv-calls", counter_type_relative);
  create("mpi-sync-time",  counter_type_relative);
  create("mpi-sync-procs", counter_type_relative);
  create("mpi-sync-calls", counter_type_relative);
}

//----------------------------------------------------------------------

CountersMpi::~CountersMpi() throw ()
{
}

//----------------------------------------------------------------------

CountersMpi::CountersMpi(const CountersMpi & counters) throw ()
  : CountersUser(counters)
/// @param     counters  Object being copied
{
}

//----------------------------------------------------------------------

CountersMpi & CountersMpi::operator= (const CountersMpi & counters) throw ()
/// @param     counters  Source object of the assignment
/// @return    The target assigned object
{
  return *this;
}

//======================================================================

long long * CountersMpi::start_()
{
  long long * counters = CountersUser::start_();

  assign_all_();

  return counters;

}

//----------------------------------------------------------------------

void CountersMpi::stop_(long long * counters)
{
  assign_all_();

  CountersUser::stop_(counters);
}

//----------------------------------------------------------------------

void CountersMpi::assign_all_()
{
  assign("mpi-time",       mpi_time_);
  assign("mpi-calls",      mpi_calls_);
  assign("mpi-send-time",  mpi_send_time_);
  assign("mpi-send-bytes", mpi_send_bytes_);
  assign("mpi-send-calls", mpi_send_calls_);
  assign("mpi-recv-time",  mpi_recv_time_);
  assign("mpi-recv-bytes", mpi_recv_bytes_);
  assign("mpi-recv-calls", mpi_recv_calls_);
  assign("mpi-sync-time",  mpi_sync_time_);
  assign("mpi-sync-procs", mpi_sync_procs_);
  assign("mpi-sync-calls", mpi_sync_calls_);
}

//----------------------------------------------------------------------

void CountersMpi::initialize()
{
  int i;

  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi_);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip_mpi_);

  long id_node = gethostid();

  int * id_node_list = new int [np_mpi_];

  for (i=0; i<np_mpi_; i++) {
    id_node_list[i] = -1;
  }

  // Get all id_node values from all processes, eg 8323328, ...
  MPI_Allgather (&id_node,1,MPI_INT,id_node_list,1,MPI_INT,MPI_COMM_WORLD);

  // Count number of id's less than this one
  int ip_node_mul = 0;
  for (i=0; i<np_mpi_; i++) {
    if (id_node > id_node_list[i]) ++ ip_node_mul;
  }

  // Get all ip_node_mul values from all processes, eg. 8, 4, 4, 0, 8, ...
  MPI_Allgather (&ip_node_mul,1,MPI_INT,id_node_list,1,MPI_INT,MPI_COMM_WORLD);

  // Count number of ip_node_mul == 0 to get np_thread_
  int np_thread = 0;
  for (i=0; i<np_mpi_; i++) {
    if (id_node_list[i] == 0) ++ np_thread;
  }

  // Compute ip_node, np_node 

  ip_node_ = ip_node_mul / np_thread;
  // Complicated, but just ceil (np_mpi / np_thread)
  np_node_ = (np_mpi_ + (np_thread - 1 - ((np_mpi_-1)%np_thread))) / np_thread;


//   // Compute ip_thread.  Note ip_mpi upper loop limit

//   ip_thread_ = 0;
//   for (i=0; i<ip_mpi_; i++) {
//     if (id_node_list[i] == ip_node_mul) ++ ip_thread_;
//   }
}

//======================================================================

int MPI_Send (void* b, int n, MPI_Datatype t, int d,int g, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Send(b,n,t,d,g,c);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);

  return result;

}

//----------------------------------------------------------------------

int MPI_Bsend(void* b, int n, MPI_Datatype t, int d,int g, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Bsend(b,n,t,d,g,c);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);

  return result;
}

//----------------------------------------------------------------------

int MPI_Rsend(void* b, int n, MPI_Datatype t, int d,int g, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Rsend(b,n,t,d,g,c);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);

  return result;
}

//----------------------------------------------------------------------

int MPI_Ssend(void* b, int n, MPI_Datatype t, int d,int g, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Ssend(b,n,t,d,g,c);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);

  return result;
}

//----------------------------------------------------------------------

int MPI_Isend (void* b, int n, MPI_Datatype t, int d, int g, MPI_Comm c, 
	       MPI_Request *r)
{
  CountersMpi::timer_start();

  int result = PMPI_Isend(b,n,t,d,g,c,r);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);


  return result;
}

//----------------------------------------------------------------------

int MPI_Ibsend(void* b, int n, MPI_Datatype t, int d, int g, MPI_Comm c, 
	       MPI_Request *r)
{
  CountersMpi::timer_start();

  int result = PMPI_Ibsend(b,n,t,d,g,c,r);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);


  return result;
}

//----------------------------------------------------------------------

int MPI_Irsend(void* b, int n, MPI_Datatype t, int d, int g, MPI_Comm c, 
	       MPI_Request *r)
{
  CountersMpi::timer_start();

  int result = PMPI_Irsend(b,n,t,d,g,c,r);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);

  return result;
}

//----------------------------------------------------------------------

int MPI_Issend(void* b, int n, MPI_Datatype t, int d, int g, MPI_Comm c, 
	       MPI_Request *r)
{
  CountersMpi::timer_start();

  int result = PMPI_Issend(b,n,t,d,g,c,r);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);

  return result;
}

//----------------------------------------------------------------------

int MPI_Irecv(void* b, int n, MPI_Datatype t, int s, int g, MPI_Comm c, 
	      MPI_Request *r)
{
  CountersMpi::timer_start();

  int result = PMPI_Irecv(b,n,t,s,g,c,r);

  CountersMpi::timer_stop();

  CountersMpi::count_recv(n,t);

  return result;
}

//----------------------------------------------------------------------

int MPI_Recv (void* b, int n, MPI_Datatype t, int s, int g, MPI_Comm c, 
	      MPI_Status *st)
{
  CountersMpi::timer_start();

  int result = PMPI_Recv(b,n,t,s,g,c,st);

  CountersMpi::timer_stop();

  CountersMpi::count_recv(n,t);

  return result;
}

//----------------------------------------------------------------------

int MPI_Sendrecv(void *bs, int ns, MPI_Datatype ts, int d, int gs, 
                 void *br, int nr, MPI_Datatype tr, int s, int gr, 
                 MPI_Comm c, MPI_Status *st)
{
  CountersMpi::timer_start();

  int result = PMPI_Sendrecv(bs,ns,ts,d,gs,br,nr,tr,s,gr,c,st);

  CountersMpi::timer_stop();

  CountersMpi::count_send(ns,ts);
  CountersMpi::count_recv(nr,tr);

  return result;
}

//----------------------------------------------------------------------

int MPI_Sendrecv_replace(void* b, int n, MPI_Datatype t,
			 int d, int gs, int s, int gr,
			 MPI_Comm c, MPI_Status *st)
{
  CountersMpi::timer_start();

  int result = PMPI_Sendrecv_replace(b,n,t,d,gs,s,gr,c,st);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);
  CountersMpi::count_recv(n,t);

  return result;
}

//----------------------------------------------------------------------

int MPI_Allgather(void* bs, int ns, MPI_Datatype ts,
		  void* br, int nr, MPI_Datatype tr,
		  MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Allgather(bs,ns,ts,br,nr,tr,c);

  CountersMpi::timer_stop();

  int np;
  MPI_Comm_size(c,&np);

  CountersMpi::count_send(ns,ts,np-1);
  CountersMpi::count_recv(nr,tr,np-1);

  return result;
}

//----------------------------------------------------------------------

int MPI_Allgatherv(void* bs, int ns, MPI_Datatype ts,
		   void* br, int *nrs, int *ds,
		   MPI_Datatype tr, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Allgatherv(bs,ns,ts,br,nrs,ds,tr,c);

  CountersMpi::timer_stop();

  int ip,np;
  MPI_Comm_size(c,&np);
  MPI_Comm_rank(c,&ip);

  CountersMpi::count_send(ns,ts,np-1);
  for (int i=0; i<np; i++) {
    if (i != ip)  CountersMpi::count_recv(nrs[i],tr);
  }

  return result;
}

//----------------------------------------------------------------------

int MPI_Allreduce(void* bs, void* br, int n,
		  MPI_Datatype t, MPI_Op op, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Allreduce(bs,br,n,t,op,c);

  CountersMpi::timer_stop();

  CountersMpi::count_send(n,t);
  CountersMpi::count_recv(n,t);
  CountersMpi::count_sync(c);

  return result;
}

//----------------------------------------------------------------------

int MPI_Alltoall(void* bs, int ns, MPI_Datatype ts,
		 void* br, int nr, MPI_Datatype tr,
		 MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Alltoall(bs,ns,ts,br,nr,tr,c);

  CountersMpi::timer_stop();

  int np;
  MPI_Comm_size(c,&np);

  CountersMpi::count_send(ns,ts,np-1);
  CountersMpi::count_recv(nr,tr,np-1);

  return result;
}

//----------------------------------------------------------------------

int MPI_Alltoallv(void* bs, int *nss, int *dss, MPI_Datatype ts, 
		  void* br, int *nrs, int *drs, MPI_Datatype tr, 
		  MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Alltoallv(bs,nss,dss,ts,br,nrs,drs,tr,c);

  CountersMpi::timer_stop();

  int ip,np;
  MPI_Comm_rank(c,&ip);
  MPI_Comm_size(c,&np);

  for (int i=0; i<np; i++) {
    if (i != ip) {
      CountersMpi::count_send(nss[i],ts);
      CountersMpi::count_recv(nrs[i],tr);
    }
  }

  return result;
}

//----------------------------------------------------------------------

int MPI_Alltoallw(void *bs, int * nss, int * dss, MPI_Datatype * tss, 
		  void *br, int * nrs, int * drs, MPI_Datatype * trs, 
		  MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Alltoallw(bs,nss,dss,tss,br,nrs,drs,trs,c);

  CountersMpi::timer_stop();

  int ip,np;
  MPI_Comm_rank(c,&ip);
  MPI_Comm_size(c,&np);

  for (int i=0; i<np; i++) {
    if (i != ip) {
      CountersMpi::count_send(nss[i],tss[i]);
      CountersMpi::count_recv(nrs[i],tss[i]);
    }
  }

  return result;
}

//----------------------------------------------------------------------

int MPI_Barrier(MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Barrier(c);

  CountersMpi::timer_stop();

  CountersMpi::count_sync(c);

  return result;
}

//----------------------------------------------------------------------

int MPI_Bcast(void* b, int n, MPI_Datatype t, int r, MPI_Comm c )
{
  CountersMpi::timer_start();

  int result = PMPI_Bcast(b,n,t,r,c);

  CountersMpi::timer_stop();

  int ip,np;
  MPI_Comm_size(c,&np);
  MPI_Comm_rank(c,&ip);
  if (ip == r) {
    CountersMpi::count_send(n,t);
  } else {
    CountersMpi::count_recv(n,t,np-1);
  }

  return result;
}

//----------------------------------------------------------------------

int MPI_Gather(void* bs, int ns, MPI_Datatype ts,
	       void* br, int nr, MPI_Datatype tr, int r,
	       MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Gather(bs,ns,ts,br,nr,tr,r,c);

  CountersMpi::timer_stop();

  int ip,np;
  MPI_Comm_size(c,&np);
  MPI_Comm_rank(c,&ip);
  if (ip == r) {
    CountersMpi::count_recv(nr,ts,np-1);
  } else {
    CountersMpi::count_send(ns,tr);
  }

  return result;
}

//----------------------------------------------------------------------

int MPI_Gatherv(void* bs, int ns, MPI_Datatype ts,
		void* br, int *nrs, int *ds,
		MPI_Datatype tr, int r, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Gatherv(bs,ns,ts,br,nrs,ds,tr,r,c);

  CountersMpi::timer_stop();

  int ip,np;
  MPI_Comm_size(c,&np);
  MPI_Comm_rank(c,&ip);
  if (ip == r) {
    for (int i=0; i<np; i++) {
      if (i != r) CountersMpi::count_recv(nrs[i],tr);
    }
  } else {
    CountersMpi::count_send(ns,ts);
  }

  return result;
}

//----------------------------------------------------------------------

int MPI_Reduce(void* bs, void* br, int n,
	       MPI_Datatype t, MPI_Op op, int r, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Reduce(bs,br,n,t,op,r,c);
 
  CountersMpi::timer_stop();

  int ip,np;
  MPI_Comm_size(c,&np);
  MPI_Comm_rank(c,&ip);

  if (ip==r) {
    CountersMpi::count_recv(n,t,np-1);
  } else {
    CountersMpi::count_send(n,t);
  }

  return result;
}

//----------------------------------------------------------------------

int MPI_Reduce_scatter(void* bs, void* br, int *nrs,
		       MPI_Datatype t, MPI_Op op, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Reduce_scatter(bs,br,nrs,t,op,c);

  CountersMpi::timer_stop();
 
  int ip,np;
  MPI_Comm_size(c,&np);
  MPI_Comm_rank(c,&ip);

  for (int i=0; i<np; i++) {
    if (i != ip) {
      CountersMpi::count_recv(nrs[i],t);
      CountersMpi::count_send(nrs[i],t);
    }
  }
  return result;
}

//----------------------------------------------------------------------

int MPI_Scatter(void* bs, int ns, MPI_Datatype ts,
		void* br, int nr, MPI_Datatype tr, int r,
		MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Scatter(bs,ns,ts,br,nr,tr,r,c);
 
  CountersMpi::timer_stop();

  int ip,np;
  MPI_Comm_size(c,&np);
  MPI_Comm_rank(c,&ip);

  if (ip==r) {
    CountersMpi::count_send(ns,ts,np-1);
  } else {
    CountersMpi::count_recv(nr,tr);
  }

  return result;
}

//----------------------------------------------------------------------

int MPI_Scatterv(void* bs, int *nss, int *ds, MPI_Datatype ts, 
		 void* br, int nr, MPI_Datatype tr, int r, MPI_Comm c)
{
  CountersMpi::timer_start();

  int result = PMPI_Scatterv(bs,nss,ds,ts,br,nr,tr,r,c);
 
  CountersMpi::timer_stop();

  int ip,np;
  MPI_Comm_size(c,&np);
  MPI_Comm_rank(c,&ip);

  if (ip == r) {
    for (int i=0; i<np; i++) {
      if (i != r) CountersMpi::count_send(nss[i],ts);
    }
  } else {
    CountersMpi::count_recv(nr,tr);
  }

  return result;
}

//----------------------------------------------------------------------

int MPI_Wait (MPI_Request * r, MPI_Status * s)
{
  CountersMpi::timer_start();
  int result = PMPI_Wait(r, s);
  CountersMpi::timer_stop();
  return result;
}

//----------------------------------------------------------------------

int MPI_Waitany (int n,   MPI_Request * r, 
		 int * i, MPI_Status * s)
{
  CountersMpi::timer_start();
  int result = PMPI_Waitany(n,r,i,s);
  CountersMpi::timer_stop();
  return result;
}

//----------------------------------------------------------------------

int MPI_Waitall (int n, MPI_Request * r, MPI_Status * s)
{
  CountersMpi::timer_start();
  int result = PMPI_Waitall(n,r, s);
  CountersMpi::timer_stop();
  return result;
}

//----------------------------------------------------------------------

int MPI_Waitsome (int  nr,  MPI_Request * r, 
		  int *ns, int *is, MPI_Status * s)
{
  CountersMpi::timer_start();
  int result = PMPI_Waitsome(nr, r, ns, is, s);
  CountersMpi::timer_stop();
  return result;
}

//======================================================================
// PERSISTENT REQUESTS NOT CURRENTLY HANDLED
//======================================================================

int MPI_Cancel(MPI_Request *request)
{
  fprintf (stderr,"WARNING: MPI_Cancel not traced\n");
  int result = MPI_Cancel(request);
  return result;
}

//----------------------------------------------------------------------

int MPI_Start(MPI_Request *request)
{
  fprintf (stderr,"WARNING: MPI_Start not traced\n");
  int result = MPI_Start(request);
  return result;
}

//----------------------------------------------------------------------

int MPI_Startall(int n, MPI_Request *array_of_requests)
{
  fprintf (stderr,"WARNING: MPI_Startall not traced\n");
  int result = MPI_Startall(n,array_of_requests);
  return result;
}

//----------------------------------------------------------------------

int MPI_Scan(void* sendb, void* br, int n, MPI_Datatype datatype, 
	     MPI_Op op, MPI_Comm comm )
{
  fprintf (stderr,"WARNING: MPI_Scan not traced\n");
  int result = MPI_Scan(sendb,br,n,datatype,op,comm);
  return result;
}

//----------------------------------------------------------------------

int MPI_Exscan(void *sendb, void *br, int n, MPI_Datatype datatype,
	       MPI_Op op, MPI_Comm comm)
{
  fprintf (stderr,"WARNING: MPI_Exscan not traced\n");
  int result = MPI_Exscan(sendb,br,n,datatype,op,comm);
  return result;
}
}
#endif

