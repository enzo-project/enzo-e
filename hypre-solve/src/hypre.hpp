//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Hypre class include file

/**
 * 
 * @file      hypre.hpp
 * @brief     Declaration of the Hypre class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
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

  double               resid_;   // Solver residual
  int                  iter_;    // Solver iterations

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


  int    iterations () { return iter_; };
  double residual () { return resid_; };

private:

  // init_graph() functions

  void init_graph_nonstencil_ (Grid & grid)
  { init_nonstencil_ (grid, "graph"); };

  // init_matrix() functions

  void init_matrix_stencil_    (Grid & grid);
  void init_matrix_clear_      (int part);
  void init_matrix_nonstencil_ (Grid & grid)
  { init_nonstencil_ (grid, "matrix"); };

  void init_nonstencil_ (Grid & grid, std::string phase);
  
  // init_vector() functions

  Scalar init_vector_points_  (Hierarchy            & hierarchy,
			       std::vector<Point *> & points);
  Scalar init_vector_spheres_ (Hierarchy             & hierarchy,
			       std::vector<Sphere *> & spheres);		
  Scalar init_vector_density_ (Hierarchy             & hierarchy);		

  // solve() functions

  void solve_fac_      (Hierarchy & hierarchy, int itmax, double restol);
  void solve_bicgstab_ (Hierarchy & hierarchy, int itmax, double restol);
  void solve_pfmg_     (Hierarchy & hierarchy, int itmax, double restol);

};
