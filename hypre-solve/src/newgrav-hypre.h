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

  Parameters           *parameters_; // Pointer to parameters
  Hierarchy            *hierarchy_;  // Pointer to the hierarchy

  double               resid_;   // Solver residual
  int                  iter_;    // Solver iterations

public:

  Hypre (Hierarchy  & hierarchy,
	 Parameters & parameters);

  ~Hypre ();

  void init_hierarchy (Mpi        & mpi);
  void init_stencil   ();
  void init_graph     ();
  void init_linear    (std::vector<Point *>  points);
  void solve          ();
  void evaluate       ();


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

  Scalar init_vector_points_  (std::vector<Point *> & points);
  Scalar init_vector_density_ (std::string             file_prefix,
			       bool is_packed);

  // solve() functions

  void solve_fac_      (int itmax, double restol);
  void solve_bicgstab_ (int itmax, double restol);
  void solve_pfmg_     (int itmax, double restol);

};
