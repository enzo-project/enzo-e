//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Mpi class header file

/**
 * 
 * @file      mpi.hpp
 * @brief     Declaration of the Mpi class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 */

class Mpi {
  
public:

  Mpi ()
    : comm_(0),
      np_(1),
      ip_(0)
  {
  }

  Mpi (int * argc, char ***argv)
    : comm_(MPI_COMM_WORLD)
  {
    MPI_Init (argc,argv);

    MPI_Comm_size (comm_, &np_);
    MPI_Comm_rank (comm_, &ip_);
  }

  ~Mpi ()
  {
    //    MPI_Finalize (); // MPI_Finalize() seems to be called at program exit
    //                        complains if included
  };
  bool is_root () throw () {return ip_ == 0;};
  int ip () throw () {return ip_;};
  int np () throw () {return np_;};
  void barrier () throw () { MPI_Barrier (comm_); };
    
private:

  MPI_Comm comm_;
  int np_;   // Number of processors
  int ip_;   // Rank of this processor
};



extern Mpi * pmpi;


