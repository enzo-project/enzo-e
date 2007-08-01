//======================================================================
//
//        File: mpi.hpp
//
//     Summary: Communication routines
//
// Description: Basic communication routines and data
//
//      Author: James Bordner <jobordner@ucsd.edu>
//
//        Date: 2007-04-10
//
//======================================================================

class Mpi {
  
public:

  Mpi (int * argc, char ***argv);
  int ip () throw () {return ip_;};
  int np () throw () {return np_;};
    
private:

  MPI_Comm comm_;
  int np_;   // Number of processors
  int ip_;   // Rank of this processor
};


