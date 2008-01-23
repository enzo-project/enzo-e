//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Hypre class include file

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

class Hypre {

private:

  HYPRE_SStructGrid    grid_;    // hypre grid
  HYPRE_SStructGraph   graph_;   // hypre graph
  HYPRE_SStructStencil stencil_; // hypre stencil
  HYPRE_SStructMatrix  A_;       // hypre matrix
  HYPRE_SStructVector  B_;       // hypre vector right-hand side
  HYPRE_SStructVector  X_;       // hypre vector solution
  HYPRE_SStructSolver  solver_;  // hypre solver

  Parameters           parameters_; 

public:

  Hypre (Parameters & parameters);

  void init_hierarchy (Parameters & parameters,
		       Hierarchy  & hierarchy, 
		       Mpi        & mpi);
  void init_stencil   (Hierarchy & hierarchy);
  void init_graph     (Hierarchy & hierarchy);
  void init_linear    (Parameters          & parameters,
		       Hierarchy           & hierarchy,
		       std::vector<Point *>  points,
		       std::vector<Sphere *> spheres);
  void solve          (Parameters & parameters,
		       Hierarchy & hierarchy);
  void evaluate       (Hierarchy & hierarchy);


private:

  // init_graph() functions

  void init_graph_nonstencil_ (Grid & grid)
  { init_nonstencil_ (grid, "graph"); };

  // init_matrix() functions

  void init_matrix_stencil_    (Grid & grid);
  void init_matrix_clear_      (Level & level);
  void init_matrix_nonstencil_ (Grid & grid)
  { init_nonstencil_ (grid, "matrix"); };

  void init_nonstencil_ (Grid & grid, std::string phase);
  
  // init_vector() functions

  Scalar init_vector_points_  (Hierarchy            & hierarchy,
			       std::vector<Point *> & points);
  Scalar init_vector_spheres_ (Hierarchy             & hierarchy,
			       std::vector<Sphere *> & spheres);		

  // solve() functions

  void solve_fac_      (Hierarchy & hierarchy, int itmax, double restol);
  void solve_bicgstab_ (Hierarchy & hierarchy, int itmax, double restol);
  void solve_pfmg_     (Hierarchy & hierarchy, int itmax, double restol);

};
