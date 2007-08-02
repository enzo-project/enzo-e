//======================================================================
//
//        File: 
//
//     Summary: 
//
// Description:
//
//      Author: James Bordner <jobordner@ucsd.edu>
//
//        Date: 
//
//======================================================================

class Hypre {

public:

  Hypre ();

  /// Initialize HYPRE's grid hierarchy
  void init_hierarchy (Hierarchy & hierarchy, Mpi & mpi);
  void init_stencil   (Hierarchy & hierarchy);
  void init_graph     (Hierarchy & hierarchy);
  void init_matrix    (Hierarchy & hierarchy);
  void init_rhs       (Hierarchy & hierarchy);
  void init_solver    (Hierarchy & hierarchy);
  void solve          (Hierarchy & hierarchy);
  void evaluate       (Hierarchy & hierarchy);


private:

  HYPRE_SStructStencil stencil_;

};
