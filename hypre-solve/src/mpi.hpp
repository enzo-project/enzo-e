//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Mpi class header file

/**
 * 
 * Basic communication routines and data
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-04-10
 *
 */

class Mpi {
  
public:

  Mpi ();
  Mpi (int * argc, char ***argv);
  ~Mpi ();
  bool is_root () throw () {return ip_ == 0;};
  int ip () throw () {return ip_;};
  int np () throw () {return np_;};
  void barrier () throw () { MPI_Barrier (comm_); };
    
private:

  MPI_Comm comm_;
  int np_;   // Number of processors
  int ip_;   // Rank of this processor
};



