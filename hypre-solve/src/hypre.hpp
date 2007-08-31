
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

  HYPRE_SStructGrid    grid_;    // Struct for hypre grid
  HYPRE_SStructGraph   graph_;   // Struct for hypre graph
  HYPRE_SStructStencil stencil_; // Struct for hypre stencil
  HYPRE_SStructMatrix  A_;       // Struct for hypre matrix
  HYPRE_SStructVector  B_;       // Struct for hypre vector right-hand side
  HYPRE_SStructVector  X_;       // Struct for hypre vector solution
  HYPRE_SStructSolver  solver_;  // Struct for hypre solver

public:

  Hypre ();

  void init_hierarchy (Parameters & parameters,
		       Hierarchy  & hierarchy, 
		       Mpi        & mpi);
  void init_stencil   (Hierarchy & hierarchy);
  void init_graph     (Hierarchy & hierarchy);
  void init_linear    (Parameters          & parameters,
		       Hierarchy           & hierarchy,
		       std::vector<Point *>  points,
		       std::vector<Sphere *> spheres);
  void solve          (Hierarchy & hierarchy);
  void evaluate       (Hierarchy & hierarchy);


private:

  // init_graph() functions

  void init_graph_children_            (Grid & grid);
  void init_graph_parent_              (Hierarchy & hierarchy,
				        Grid & grid);
  void init_graph_neighbors_children_  (Grid & grid);
  void init_graph_parents_neighbor_    (Grid & grid);

  // init_matrix() functions

  void init_matrix_stencil_            (Grid & grid);
  void init_matrix_clear_              (Grid & grid);
  void init_matrix_children_           (Grid & grid);
  void init_matrix_parent_             (Hierarchy & hierarchy,
                                        Grid & grid);
  void init_matrix_neighbors_children_ (Grid & grid);
  void init_matrix_parents_neighbor_   (Grid & grid);

  // init_vector() functions

  Scalar init_vector_points_             (Hierarchy            & hierarchy,
					  std::vector<Point *> & points);
  Scalar init_vector_spheres_            (Hierarchy             & hierarchy,
					  std::vector<Sphere *> & spheres);		

  // solve() functions

  void solve_fac_                      (Hierarchy & hierarchy);
  void solve_pfmg_                     (Hierarchy & hierarchy);

};
